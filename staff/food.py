"""
FOOD ATLAS — Score every food by whose side it's on in the war.

For each food compound: is it antifungal (Red Queen), profungal (Lady), or neutral?
For each food: sum the war scores of all its compounds.

Data sources:
  - FooDB (28,000 compounds, 1,000 foods): foodb.ca/downloads
  - Published antifungal MIC/IC50 data
  - Ergosterol pathway inhibitors
  - Prostaglandin synthesis data

Scoring:
  POSITIVE = helps the host (antifungal, immune-boosting, thermal wall supporting)
  NEGATIVE = helps the Lady (sugar, immune-suppressing, wall-degrading)
  ZERO = neutral

Usage:
    from staff.food import FoodAtlas

    atlas = FoodAtlas()
    atlas.score_food("garlic")       # → high positive (allicin)
    atlas.score_food("white sugar")  # → negative (pure carbon for her)
    atlas.score_food("oregano")      # → high positive (carvacrol + thymol)
    atlas.compare("garlic", "bread") # → which one fights her better
"""

import json
from pathlib import Path
from typing import Dict, Optional, List, Tuple

# ── KNOWN ANTIFUNGAL COMPOUNDS (published data) ────────────────────

ANTIFUNGAL_COMPOUNDS = {
    # Compound: {mechanism, targets, potency (1-10), sources}
    'allicin': {
        'mechanism': 'Membrane disruption. Kills 20+ Candida + Aspergillus species.',
        'targets': ['Candida', 'Aspergillus', 'Cryptococcus', 'Saccharomyces'],
        'potency': 9,
        'food_sources': ['garlic', 'onion', 'leek', 'shallot', 'chives'],
        'citation': 'Frontiers Microbiol 2022; PubMed 15743425',
    },
    'carvacrol': {
        'mechanism': 'Mimics calcium stress. Inhibits TOR pathway. Membrane disruption.',
        'targets': ['Candida', 'Aspergillus', 'broad spectrum'],
        'potency': 8,
        'food_sources': ['oregano', 'thyme', 'savory', 'marjoram'],
        'citation': 'PMC2981246',
    },
    'thymol': {
        'mechanism': 'Membrane disruption. Synergistic with carvacrol.',
        'targets': ['Candida', 'Aspergillus', 'dermatophytes'],
        'potency': 7,
        'food_sources': ['thyme', 'oregano', 'basil'],
        'citation': 'PMC2981246',
    },
    'caprylic_acid': {
        'mechanism': 'Integrates into lipid bilayer. Increases permeability. Cell lysis.',
        'targets': ['Candida'],
        'potency': 7,
        'food_sources': ['coconut oil', 'palm kernel oil', 'butter', 'human breast milk'],
        'citation': 'Multiple Candida studies',
    },
    'capric_acid': {
        'mechanism': 'Same membrane disruption as caprylic. Medium-chain fatty acid.',
        'targets': ['Candida'],
        'potency': 6,
        'food_sources': ['coconut oil', 'palm kernel oil', 'goat milk'],
        'citation': 'PMID 12766136',
    },
    'lauric_acid': {
        'mechanism': 'Membrane disruption. Most abundant MCFA in coconut oil.',
        'targets': ['Candida', 'broad spectrum'],
        'potency': 7,
        'food_sources': ['coconut oil', 'palm kernel oil', 'human breast milk'],
        'citation': 'Multiple studies',
    },
    'curcumin': {
        'mechanism': 'Disrupts ergosterol synthesis. Membrane damage. ROS generation.',
        'targets': ['Candida', 'Aspergillus', 'Penicillium'],
        'potency': 6,
        'food_sources': ['turmeric', 'curry powder', 'mustard'],
        'citation': 'Multiple antifungal studies',
    },
    'cinnamaldehyde': {
        'mechanism': 'Membrane disruption. Inhibits cell wall synthesis.',
        'targets': ['Candida', 'Aspergillus'],
        'potency': 7,
        'food_sources': ['cinnamon', 'cassia'],
        'citation': 'Multiple studies',
    },
    'eugenol': {
        'mechanism': 'Membrane disruption. Ergosterol binding.',
        'targets': ['Candida', 'Aspergillus', 'Penicillium'],
        'potency': 6,
        'food_sources': ['cloves', 'allspice', 'basil', 'cinnamon', 'nutmeg'],
        'citation': 'PMID 15928026',
    },
    'acetic_acid': {
        'mechanism': 'pH disruption. Membrane permeabilization at low pH.',
        'targets': ['Candida', 'broad spectrum'],
        'potency': 4,
        'food_sources': ['vinegar', 'fermented vegetables', 'kombucha'],
        'citation': 'Multiple studies',
    },
    'tannins': {
        'mechanism': 'Protein precipitation. Enzyme inhibition. Iron chelation.',
        'targets': ['Candida', 'broad spectrum'],
        'potency': 5,
        'food_sources': ['black tea', 'green tea', 'red wine', 'pomegranate',
                         'cranberry', 'walnut', 'dark chocolate'],
        'citation': 'Multiple studies',
    },
    'lactoferrin': {
        'mechanism': 'Iron sequestration. Starves fungi of essential iron.',
        'targets': ['Candida', 'Aspergillus', 'broad spectrum'],
        'potency': 7,
        'food_sources': ['breast milk', 'colostrum', 'raw milk', 'whey'],
        'citation': 'PMID 15558068',
    },
    'ajoene': {
        'mechanism': 'Thiosulfinate. Membrane disruption. Biofilm inhibition.',
        'targets': ['Candida', 'Aspergillus'],
        'potency': 8,
        'food_sources': ['garlic (aged/crushed)'],
        'citation': 'Nature Sci Rep 2018',
    },
    'berberine': {
        'mechanism': 'AMPK activator. Membrane disruption. DNA intercalation.',
        'targets': ['Candida', 'Cryptococcus'],
        'potency': 7,
        'food_sources': ['goldenseal', 'Oregon grape', 'barberry', 'turmeric tree'],
        'citation': 'Multiple studies',
    },
    'propolis_compounds': {
        'mechanism': 'Mixed flavonoids + phenolics. Membrane disruption.',
        'targets': ['Candida', 'Aspergillus', 'dermatophytes'],
        'potency': 6,
        'food_sources': ['bee propolis', 'raw honey'],
        'citation': 'Multiple studies',
    },
    'beta_glucan_dietary': {
        'mechanism': 'Immune stimulation. Activates Dectin-1 (anti-fungal receptor).',
        'targets': ['Immune boost against all fungi'],
        'potency': 5,
        'food_sources': ['oats', 'barley', 'mushrooms', 'nutritional yeast'],
        'citation': 'Multiple immunology studies',
    },
    'sulforaphane': {
        'mechanism': 'NRF2 activation. Antioxidant. Some antifungal activity.',
        'targets': ['Candida', 'broad spectrum'],
        'potency': 4,
        'food_sources': ['broccoli', 'broccoli sprouts', 'cauliflower', 'kale',
                         'brussels sprouts', 'cabbage'],
        'citation': 'Multiple studies',
    },
}

# ── PROFUNGAL COMPOUNDS (help the Lady) ────────────────────────────

PROFUNGAL_COMPOUNDS = {
    'glucose': {
        'mechanism': 'Primary carbon source. She eats it directly.',
        'potency': -8,
        'food_sources': ['white sugar', 'candy', 'soft drinks', 'fruit juice',
                         'white bread', 'white rice', 'corn syrup', 'honey'],
    },
    'fructose': {
        'mechanism': 'Carbon source. Fermentable sugar. Feeds Candida.',
        'potency': -7,
        'food_sources': ['fruit juice', 'high-fructose corn syrup', 'agave',
                         'soft drinks', 'dried fruit'],
    },
    'refined_starch': {
        'mechanism': 'Rapidly converted to glucose. Feeds her.',
        'potency': -5,
        'food_sources': ['white bread', 'white rice', 'pasta', 'crackers',
                         'processed cereals', 'pastries'],
    },
    'alcohol_ethanol': {
        'mechanism': 'Immunosuppressive. Disrupts gut barrier. Feeds yeast directly.',
        'potency': -7,
        'food_sources': ['beer', 'wine', 'spirits', 'cocktails'],
    },
    'aflatoxin': {
        'mechanism': 'HER chemical weapon. Produced by Aspergillus on stored crops. IARC Group 1 carcinogen.',
        'potency': -10,
        'food_sources': ['contaminated peanuts', 'contaminated corn',
                         'contaminated tree nuts', 'contaminated grains'],
    },
}


class FoodAtlas:
    """Score foods by whose side they're on in the war."""

    def __init__(self):
        self.antifungal = ANTIFUNGAL_COMPOUNDS
        self.profungal = PROFUNGAL_COMPOUNDS
        self._food_index = self._build_food_index()

    def _build_food_index(self) -> Dict[str, Dict]:
        """Build reverse index: food → list of compounds."""
        index = {}
        for name, data in self.antifungal.items():
            for food in data['food_sources']:
                food_lower = food.lower()
                if food_lower not in index:
                    index[food_lower] = {'antifungal': [], 'profungal': [], 'score': 0}
                index[food_lower]['antifungal'].append((name, data['potency']))
                index[food_lower]['score'] += data['potency']

        for name, data in self.profungal.items():
            for food in data['food_sources']:
                food_lower = food.lower()
                if food_lower not in index:
                    index[food_lower] = {'antifungal': [], 'profungal': [], 'score': 0}
                index[food_lower]['profungal'].append((name, data['potency']))
                index[food_lower]['score'] += data['potency']

        return index

    def score_food(self, food_name: str) -> Dict:
        """Score a food: positive = fights her, negative = feeds her."""
        food_lower = food_name.lower()

        # Exact match
        if food_lower in self._food_index:
            return self._food_index[food_lower]

        # Partial match
        matches = [k for k in self._food_index if food_lower in k or k in food_lower]
        if matches:
            best = max(matches, key=lambda k: abs(self._food_index[k]['score']))
            result = self._food_index[best]
            result['matched'] = best
            return result

        return {'antifungal': [], 'profungal': [], 'score': 0, 'unknown': True}

    def compare(self, food_a: str, food_b: str) -> str:
        """Compare two foods: which fights her better?"""
        a = self.score_food(food_a)
        b = self.score_food(food_b)
        diff = a['score'] - b['score']
        if diff > 0:
            return f"{food_a} (score {a['score']}) fights her better than {food_b} (score {b['score']})"
        elif diff < 0:
            return f"{food_b} (score {b['score']}) fights her better than {food_a} (score {a['score']})"
        else:
            return f"Equal: both score {a['score']}"

    def top_antifungal(self, n: int = 20) -> List[Tuple[str, int]]:
        """Return top n foods that fight her."""
        ranked = sorted(self._food_index.items(), key=lambda x: -x[1]['score'])
        return [(name, data['score']) for name, data in ranked[:n] if data['score'] > 0]

    def top_profungal(self, n: int = 20) -> List[Tuple[str, int]]:
        """Return top n foods that feed her."""
        ranked = sorted(self._food_index.items(), key=lambda x: x[1]['score'])
        return [(name, data['score']) for name, data in ranked[:n] if data['score'] < 0]

    def war_report(self):
        """Print the full food war report."""
        print("=" * 60)
        print("  FOOD WAR REPORT — Whose Side Is Your Diet On?")
        print("=" * 60)

        print("\n  TOP ANTIFUNGAL FOODS (fight the Lady):")
        for name, score in self.top_antifungal(15):
            bar = '+' * score
            print(f"    {score:+3d}  {name:30s}  {bar}")

        print("\n  TOP PROFUNGAL FOODS (feed the Lady):")
        for name, score in self.top_profungal(15):
            bar = '-' * abs(score)
            print(f"    {score:+3d}  {name:30s}  {bar}")

        print(f"\n  Total indexed foods: {len(self._food_index)}")
        antifungal_count = sum(1 for d in self._food_index.values() if d['score'] > 0)
        profungal_count = sum(1 for d in self._food_index.values() if d['score'] < 0)
        print(f"  Antifungal: {antifungal_count}")
        print(f"  Profungal: {profungal_count}")


def main():
    atlas = FoodAtlas()
    atlas.war_report()

    import sys
    if len(sys.argv) > 1:
        food = ' '.join(sys.argv[1:])
        print(f"\n  Scoring: {food}")
        result = atlas.score_food(food)
        print(f"  Score: {result['score']:+d}")
        if result['antifungal']:
            print(f"  Antifungal compounds: {result['antifungal']}")
        if result['profungal']:
            print(f"  Profungal compounds: {result['profungal']}")


if __name__ == '__main__':
    main()
