"""
Shared BAM streaming — pure binary parser, no external libraries.
Reads BGZF blocks, parses BAM records, yields alignments.
Used by all STAFF BAM tools.
"""
import struct
import zlib

SEQ_DECODE = '=ACMGRSVTWYHKDBN'

BASE_TO_BITS = {
    'A': 0b00, 'T': 0b01, 'C': 0b10, 'G': 0b11,
    'a': 0b00, 't': 0b01, 'c': 0b10, 'g': 0b11,
    'N': 0b00, 'n': 0b00,
}
BITS_TO_BASE = {0b00: 'A', 0b01: 'T', 0b10: 'C', 0b11: 'G'}


def seq_to_binary(seq):
    """Convert ATCG string to list of 2-bit values."""
    import numpy as np
    return np.array([BASE_TO_BITS.get(b, 0) for b in seq], dtype=np.uint8)


def decode_bam_seq(encoded_bytes, seq_len):
    """Decode BAM 4-bit packed sequence to ATCG string."""
    seq = []
    for i in range(seq_len):
        byte_idx = i >> 1
        if byte_idx >= len(encoded_bytes):
            break
        if i & 1:
            val = encoded_bytes[byte_idx] & 0x0F
        else:
            val = (encoded_bytes[byte_idx] >> 4) & 0x0F
        seq.append(SEQ_DECODE[val])
    return ''.join(seq)


def read_bgzf_block(fh):
    """Read one BGZF block, return decompressed bytes or None."""
    header = fh.read(18)
    if len(header) < 18:
        return None
    if header[:4] != b'\x1f\x8b\x08\x04':
        return None
    bsize = struct.unpack('<H', header[16:18])[0]
    remaining = bsize - 18 + 1
    cdata = fh.read(remaining)
    if len(cdata) < remaining:
        return None
    return zlib.decompress(cdata[:-8], -15)


def stream_bam(bam_path, min_mapq=10):
    """Stream BAM yielding (chrom, pos, mapq, cigar_ops, seq, flag) tuples."""
    with open(bam_path, 'rb') as fh:
        first_block = read_bgzf_block(fh)
        if first_block is None or first_block[:4] != b'BAM\x01':
            raise ValueError(f'Not a valid BAM: {bam_path}')

        data = first_block
        pos = 4
        l_text = struct.unpack('<i', data[pos:pos+4])[0]
        pos += 4 + l_text
        n_ref = struct.unpack('<i', data[pos:pos+4])[0]
        pos += 4

        ref_names = []
        for _ in range(n_ref):
            l_name = struct.unpack('<i', data[pos:pos+4])[0]
            pos += 4
            name = data[pos:pos+l_name-1].decode()
            pos += l_name
            l_ref = struct.unpack('<i', data[pos:pos+4])[0]
            pos += 4
            ref_names.append(name)

        buffer = data[pos:]
        n_yielded = 0

        while True:
            while len(buffer) < 4:
                block = read_bgzf_block(fh)
                if block is None:
                    return
                buffer += block

            block_size = struct.unpack('<i', buffer[:4])[0]

            while len(buffer) < block_size + 4:
                block = read_bgzf_block(fh)
                if block is None:
                    return
                buffer += block

            record = buffer[4:block_size+4]
            buffer = buffer[block_size+4:]

            if len(record) < 32:
                continue

            ref_id = struct.unpack('<i', record[0:4])[0]
            aln_pos = struct.unpack('<i', record[4:8])[0]
            l_read_name = record[8]
            mapq = record[9]
            n_cigar_op = struct.unpack('<H', record[12:14])[0]
            flag = struct.unpack('<H', record[14:16])[0]
            l_seq = struct.unpack('<i', record[16:20])[0]

            if ref_id < 0 or flag & 4 or mapq < min_mapq:
                continue

            cigar_start = 32 + l_read_name
            cigar_end = cigar_start + n_cigar_op * 4
            cigar_ops = []
            for i in range(n_cigar_op):
                val = struct.unpack('<I', record[cigar_start + i*4:cigar_start + (i+1)*4])[0]
                cigar_ops.append((val & 0xF, val >> 4))

            seq_start = cigar_end
            seq_bytes = (l_seq + 1) // 2
            seq_end = seq_start + seq_bytes
            if seq_end > len(record):
                continue

            seq = decode_bam_seq(record[seq_start:seq_end], l_seq)
            chrom = ref_names[ref_id] if ref_id < len(ref_names) else str(ref_id)

            yield chrom, aln_pos, mapq, cigar_ops, seq, flag

            n_yielded += 1
            if n_yielded % 200000 == 0:
                print(f'  {n_yielded:,} reads processed...')
