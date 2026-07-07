import re
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

st.set_page_config(page_title="Experimental Helper", page_icon="🧪", layout="centered")

# ---------------------------------------------------------------------------
# 원자량 (g/mol)
# ---------------------------------------------------------------------------
ATOMIC_WEIGHTS = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007,
    "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982,
    "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098,
    "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723,
    "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798, "Rb": 85.468,
    "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95, "Tc": 98,
    "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82,
    "Sn": 118.71, "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91,
    "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Pm": 145,
    "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.50, "Ho": 164.93,
    "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97, "Hf": 178.49, "Ta": 180.95,
    "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97,
    "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209, "At": 210, "Rn": 222,
    "Fr": 223, "Ra": 226, "Ac": 227, "Th": 232.04, "Pa": 231.04, "U": 238.03, "Np": 237,
    "Pu": 244, "Am": 243, "Cm": 247, "Bk": 247, "Cf": 251, "Es": 252, "Fm": 257, "Md": 258,
    "No": 259, "Lr": 262,
}

# ---------------------------------------------------------------------------
# 자주 쓰는 시약: (표시용 화학식, 한글 이름, RDKit SMILES)
# 화학식 파싱/몰질량 계산은 이 목록과 무관하게 어떤 유효한 화학식이든 가능하지만,
# 실제 골격 구조식은 SMILES가 등록된 화합물에 대해서만 그려진다.
# ---------------------------------------------------------------------------
REAGENTS = [
    ("NaCl", "염화나트륨", "[Na+].[Cl-]"),
    ("NaOH", "수산화나트륨", "[Na+].[OH-]"),
    ("KOH", "수산화칼륨", "[K+].[OH-]"),
    ("HCl", "염산", "Cl"),
    ("H2SO4", "황산", "OS(=O)(=O)O"),
    ("HNO3", "질산", "O[N+](=O)[O-]"),
    ("H3PO4", "인산", "OP(=O)(O)O"),
    ("NaHCO3", "탄산수소나트륨", "[Na+].OC([O-])=O"),
    ("Na2CO3", "탄산나트륨", "[Na+].[Na+].[O-]C([O-])=O"),
    ("K2CO3", "탄산칼륨", "[K+].[K+].[O-]C([O-])=O"),
    ("CaCl2", "염화칼슘", "[Ca+2].[Cl-].[Cl-]"),
    ("KCl", "염화칼륨", "[K+].[Cl-]"),
    ("NH4Cl", "염화암모늄", "[NH4+].[Cl-]"),
    ("NaNO3", "질산나트륨", "[Na+].[O-][N+](=O)[O-]"),
    ("KNO3", "질산칼륨", "[K+].[O-][N+](=O)[O-]"),
    ("AgNO3", "질산은", "[Ag+].[O-][N+](=O)[O-]"),
    ("FeCl3", "염화철(III)", "[Fe+3].[Cl-].[Cl-].[Cl-]"),
    ("FeCl2", "염화철(II)", "[Fe+2].[Cl-].[Cl-]"),
    ("CuSO4", "황산구리(II)", "[Cu+2].[O-]S([O-])(=O)=O"),
    ("CuSO4·5H2O", "황산구리(II) 오수화물", "[Cu+2].[O-]S([O-])(=O)=O.O.O.O.O.O"),
    ("MgSO4·7H2O", "황산마그네슘 칠수화물", "[Mg+2].[O-]S([O-])(=O)=O.O.O.O.O.O.O.O"),
    ("Na2SO4", "황산나트륨", "[Na+].[Na+].[O-]S([O-])(=O)=O"),
    ("Ca(OH)2", "수산화칼슘", "[Ca+2].[OH-].[OH-]"),
    ("BaCl2", "염화바륨", "[Ba+2].[Cl-].[Cl-]"),
    ("ZnSO4", "황산아연", "[Zn+2].[O-]S([O-])(=O)=O"),
    ("Al2(SO4)3", "황산알루미늄", "[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O"),
    ("Pb(NO3)2", "질산납", "[Pb+2].[O-][N+](=O)[O-].[O-][N+](=O)[O-]"),
    ("NH4NO3", "질산암모늄", "[NH4+].[O-][N+](=O)[O-]"),
    ("KMnO4", "과망간산칼륨", "[K+].[O-][Mn](=O)(=O)=O"),
    ("CaCO3", "탄산칼슘", "[Ca+2].[O-]C([O-])=O"),
    ("BaSO4", "황산바륨", "[Ba+2].[O-]S([O-])(=O)=O"),
    ("AgCl", "염화은", "[Ag+].[Cl-]"),
    ("CuCl2", "염화구리(II)", "[Cu+2].[Cl-].[Cl-]"),
    ("MgCl2", "염화마그네슘", "[Mg+2].[Cl-].[Cl-]"),
    ("MgO", "산화마그네슘", "[Mg+2].[O-2]"),
    ("CaO", "산화칼슘", "[Ca+2].[O-2]"),
    ("NaF", "불화나트륨", "[Na+].[F-]"),
    ("KBr", "브롬화칼륨", "[K+].[Br-]"),
    ("NaBr", "브롬화나트륨", "[Na+].[Br-]"),
    ("KI", "요오드화칼륨", "[K+].[I-]"),
    ("NaI", "요오드화나트륨", "[Na+].[I-]"),
    ("ZnCl2", "염화아연", "[Zn+2].[Cl-].[Cl-]"),
    ("AlCl3", "염화알루미늄", "[Al+3].[Cl-].[Cl-].[Cl-]"),
    ("H2O2", "과산화수소", "OO"),
    ("NaOCl", "차아염소산나트륨", "[Na+].[O-]Cl"),
    ("H3BO3", "붕산", "OB(O)O"),
    ("H2O", "물", "O"),
    ("CO2", "이산화탄소", "O=C=O"),
    ("CH4", "메탄", "C"),
    ("NH3", "암모니아", "N"),
    ("N2", "질소", "N#N"),
    ("O2", "산소", "O=O"),
    ("H2", "수소", "[HH]"),
    ("Cl2", "염소", "ClCl"),
    ("SO2", "이산화황", "O=S=O"),
    ("SO3", "삼산화황", "O=S(=O)=O"),
    ("H2S", "황화수소", "S"),
    ("PCl3", "삼염화인", "ClP(Cl)Cl"),
    ("BF3", "삼불화붕소", "FB(F)F"),
    ("CCl4", "사염화탄소", "ClC(Cl)(Cl)Cl"),
    ("CO", "일산화탄소", "[C-]#[O+]"),
    ("CH3COOH", "아세트산", "CC(=O)O"),
    ("C2H5OH", "에탄올", "CCO"),
    ("CH3OH", "메탄올", "CO"),
    ("C3H8O", "이소프로필알코올", "CC(C)O"),
    ("C3H6O", "아세톤", "CC(C)=O"),
    ("C6H12O6", "포도당", "OCC1OC(O)C(O)C(O)C1O"),
    ("C3H8O3", "글리세린", "OCC(O)CO"),
    ("C3H8O2", "프로필렌글리콜", "CC(O)CO"),
    ("C2H6O2", "에틸렌글리콜", "OCCO"),
    ("C6H8O7", "구연산", "OC(=O)CC(O)(CC(=O)O)C(=O)O"),
    ("CH4N2O", "요소", "NC(=O)N"),
    ("C3H6O3", "젖산", "CC(O)C(=O)O"),
    ("C2H4O3", "글리콜산", "OCC(=O)O"),
    ("C7H6O3", "살리실산", "OC(=O)c1ccccc1O"),
    ("C7H6O2", "벤조산", "OC(=O)c1ccccc1"),
    ("C7H5NaO2", "벤조산나트륨", "[Na+].[O-]C(=O)c1ccccc1"),
    ("C6H6N2O", "나이아신아마이드", "NC(=O)c1cccnc1"),
    ("CH2O", "포름알데히드", "C=O"),
    ("C7H8O", "벤질알코올", "OCc1ccccc1"),
    ("C6H6O", "페놀", "Oc1ccccc1"),
    ("C8H10O2", "페녹시에탄올", "OCCOc1ccccc1"),
    ("C8H8O3", "메틸파라벤", "COC(=O)c1ccc(O)cc1"),
    ("C10H12O3", "프로필파라벤", "CCCOC(=O)c1ccc(O)cc1"),
    ("C10H16N2O8", "EDTA", "OC(=O)CN(CC(=O)O)CCN(CC(=O)O)CC(=O)O"),
    ("C12H25NaO4S", "라우릴황산나트륨(SLS)", "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]"),
    ("C6H15NO3", "트리에탄올아민", "OCCN(CCO)CCO"),
    ("C8H10N4O2", "카페인", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"),
]

COMMON_REAGENTS = [(f, n) for f, n, _ in REAGENTS]

# 잘 알려진 시약의 실제 색을 저울 더미에 반영하기 위한 매칭 규칙 (구체적인 것을 먼저 검사)
REAGENT_COLOR_RULES = [
    ("CUSO4", "#8ec4f5", "#1c5fa8"),
    ("KMNO4", "#9b6bc9", "#4a2570"),
    ("FECL3", "#e0a659", "#8a5a1e"),
    ("AGNO3", "#f5f5f0", "#d8d8cc"),
    ("NACL", "#f7f7f2", "#e2e2d8"),
    ("NAOH", "#f7f7f2", "#e2e2d8"),
    ("KOH", "#f7f7f2", "#e2e2d8"),
    ("KCL", "#f7f7f2", "#e2e2d8"),
    ("NAHCO3", "#f7f7f2", "#e2e2d8"),
    ("NA2CO3", "#f7f7f2", "#e2e2d8"),
    ("K2CO3", "#f7f7f2", "#e2e2d8"),
    ("NH4CL", "#f7f7f2", "#e2e2d8"),
    ("NANO3", "#f7f7f2", "#e2e2d8"),
    ("KNO3", "#f7f7f2", "#e2e2d8"),
    ("MGSO4", "#f7f7f2", "#e2e2d8"),
    ("NA2SO4", "#f7f7f2", "#e2e2d8"),
    ("CA(OH)2", "#f7f7f2", "#e2e2d8"),
    ("C6H12O6", "#fbfaf5", "#e8e3d0"),
]
DEFAULT_PILE_COLORS = ("#fff7e3", "#e8cd8a")


# ---------------------------------------------------------------------------
# 화학식 파싱 (기존 클라이언트 JS 파서를 그대로 이식)
# ---------------------------------------------------------------------------
def normalize_key(formula: str) -> str:
    stripped = re.sub(r"\s+", "", formula)
    stripped = re.sub(r"[.*]", "·", stripped)
    return stripped.upper()


def _parse_group(formula: str) -> dict:
    pos = 0
    length = len(formula)

    def read_number():
        nonlocal pos
        start = pos
        while pos < length and formula[pos].isdigit():
            pos += 1
        if pos == start:
            return 1
        return int(formula[start:pos])

    def parse_expr():
        nonlocal pos
        counts = {}
        while pos < length:
            ch = formula[pos]
            if ch in "([":
                pos += 1
                inner = parse_expr()
                close = ")" if ch == "(" else "]"
                if pos >= length or formula[pos] != close:
                    raise ValueError("mismatched brackets")
                pos += 1
                n = read_number()
                for el, c in inner.items():
                    counts[el] = counts.get(el, 0) + c * n
            elif ch in ")]":
                break
            elif ch.isupper() and ch.isalpha():
                el = ch
                pos += 1
                if pos < length and formula[pos].islower():
                    el += formula[pos]
                    pos += 1
                if el not in ATOMIC_WEIGHTS:
                    raise ValueError("unknown element: " + el)
                n = read_number()
                counts[el] = counts.get(el, 0) + n
            else:
                raise ValueError("unexpected character: " + ch)
        return counts

    result = parse_expr()
    if pos != length:
        raise ValueError("unexpected character: " + formula[pos])
    return result


def parse_formula(formula: str) -> dict:
    parts = [p.strip() for p in re.split(r"[·.*]", formula) if p.strip()]
    if not parts:
        raise ValueError("empty formula")

    total = {}
    for part in parts:
        m = re.match(r"^(\d+)?(.+)$", part)
        if not m:
            raise ValueError("invalid segment: " + part)
        multiplier = int(m.group(1)) if m.group(1) else 1
        counts = _parse_group(m.group(2))
        for el, c in counts.items():
            total[el] = total.get(el, 0) + c * multiplier
    return total


def molar_mass(formula: str) -> float:
    counts = parse_formula(formula)
    return sum(ATOMIC_WEIGHTS[el] * n for el, n in counts.items())


def get_pile_colors(formula_text: str):
    key = normalize_key(formula_text)
    for keyword, top, bottom in REAGENT_COLOR_RULES:
        if keyword in key:
            return top, bottom
    return DEFAULT_PILE_COLORS


def format_number(n: float) -> str:
    if n != n or n in (float("inf"), float("-inf")):
        return "0"
    s = f"{n:,.4f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s


# ---------------------------------------------------------------------------
# 단순 이온성 화합물(양이온 1종 + 음이온 1종)은 RDKit이 서로 떨어진 원자로
# 그리면 그냥 "Na+ Cl-" 같은 텍스트로만 보이므로, 알약 모양 이온 카드로 따로 그린다.
# (cation_symbol, cation_count, cation_charge, anion_symbol, anion_count, anion_charge)
# ---------------------------------------------------------------------------
IONIC_PAIRS = {
    normalize_key("NaOH"): ("Na", 1, 1, "OH", 1, -1),
    normalize_key("KOH"): ("K", 1, 1, "OH", 1, -1),
    normalize_key("NaHCO3"): ("Na", 1, 1, "HCO3", 1, -1),
    normalize_key("Na2CO3"): ("Na", 2, 1, "CO3", 1, -2),
    normalize_key("K2CO3"): ("K", 2, 1, "CO3", 1, -2),
    normalize_key("CaCl2"): ("Ca", 1, 2, "Cl", 2, -1),
    normalize_key("NaCl"): ("Na", 1, 1, "Cl", 1, -1),
    normalize_key("KCl"): ("K", 1, 1, "Cl", 1, -1),
    normalize_key("NH4Cl"): ("NH4", 1, 1, "Cl", 1, -1),
    normalize_key("NaNO3"): ("Na", 1, 1, "NO3", 1, -1),
    normalize_key("KNO3"): ("K", 1, 1, "NO3", 1, -1),
    normalize_key("AgNO3"): ("Ag", 1, 1, "NO3", 1, -1),
    normalize_key("FeCl3"): ("Fe", 1, 3, "Cl", 3, -1),
    normalize_key("FeCl2"): ("Fe", 1, 2, "Cl", 2, -1),
    normalize_key("CuSO4"): ("Cu", 1, 2, "SO4", 1, -2),
    normalize_key("CuSO4·5H2O"): ("Cu", 1, 2, "SO4", 1, -2),
    normalize_key("MgSO4·7H2O"): ("Mg", 1, 2, "SO4", 1, -2),
    normalize_key("Na2SO4"): ("Na", 2, 1, "SO4", 1, -2),
    normalize_key("Ca(OH)2"): ("Ca", 1, 2, "OH", 2, -1),
    normalize_key("BaCl2"): ("Ba", 1, 2, "Cl", 2, -1),
    normalize_key("ZnSO4"): ("Zn", 1, 2, "SO4", 1, -2),
    normalize_key("Al2(SO4)3"): ("Al", 2, 3, "SO4", 3, -2),
    normalize_key("Pb(NO3)2"): ("Pb", 1, 2, "NO3", 2, -1),
    normalize_key("NH4NO3"): ("NH4", 1, 1, "NO3", 1, -1),
    normalize_key("KMnO4"): ("K", 1, 1, "MnO4", 1, -1),
    normalize_key("CaCO3"): ("Ca", 1, 2, "CO3", 1, -2),
    normalize_key("BaSO4"): ("Ba", 1, 2, "SO4", 1, -2),
    normalize_key("AgCl"): ("Ag", 1, 1, "Cl", 1, -1),
    normalize_key("CuCl2"): ("Cu", 1, 2, "Cl", 2, -1),
    normalize_key("MgCl2"): ("Mg", 1, 2, "Cl", 2, -1),
    normalize_key("MgO"): ("Mg", 1, 2, "O", 1, -2),
    normalize_key("CaO"): ("Ca", 1, 2, "O", 1, -2),
    normalize_key("NaF"): ("Na", 1, 1, "F", 1, -1),
    normalize_key("KBr"): ("K", 1, 1, "Br", 1, -1),
    normalize_key("NaBr"): ("Na", 1, 1, "Br", 1, -1),
    normalize_key("KI"): ("K", 1, 1, "I", 1, -1),
    normalize_key("NaI"): ("Na", 1, 1, "I", 1, -1),
    normalize_key("ZnCl2"): ("Zn", 1, 2, "Cl", 2, -1),
    normalize_key("AlCl3"): ("Al", 1, 3, "Cl", 3, -1),
    normalize_key("NaOCl"): ("Na", 1, 1, "OCl", 1, -1),
}


def _ion_pill_width(label: str) -> float:
    return max(34, len(label) * 13 + 18)


def _ion_pill_svg(cx: float, cy: float, symbol: str, count: int, charge: int) -> str:
    w = _ion_pill_width(symbol)
    h = 28
    x = cx - w / 2
    y = cy - h / 2
    charge_text = ""
    if charge:
        mag = abs(charge)
        charge_text = (str(mag) if mag > 1 else "") + ("+" if charge > 0 else "−")
    count_sup = f'<tspan font-size="10" dy="4">{count}</tspan>' if count > 1 else ""
    svg = (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{h / 2}" '
        f'fill="#e4ebff" stroke="#8ea3e6" stroke-width="1.3"/>'
        f'<text x="{cx}" y="{cy + 1}" text-anchor="middle" dominant-baseline="central" '
        f'font-size="15" font-weight="700" font-family="Arial, sans-serif" fill="#1c2942">'
        f'{symbol}{count_sup}</text>'
    )
    if charge_text:
        svg += (
            f'<text x="{x + w - 2}" y="{y - 4}" text-anchor="middle" font-size="12" '
            f'font-weight="700" font-family="Arial, sans-serif" fill="#3a5fad">{charge_text}</text>'
        )
    return svg


def build_ionic_svg(cation_symbol, cation_count, cation_charge, anion_symbol, anion_count, anion_charge) -> str:
    left_w = _ion_pill_width(cation_symbol)
    right_w = _ion_pill_width(anion_symbol)
    gap = 26
    total_w = left_w + right_w + gap
    start_x = (220 - total_w) / 2
    left_cx = start_x + left_w / 2
    right_cx = start_x + left_w + gap + right_w / 2
    cy = 80
    svg = '<svg viewBox="0 0 220 150" xmlns="http://www.w3.org/2000/svg">'
    svg += _ion_pill_svg(left_cx, cy, cation_symbol, cation_count, cation_charge)
    svg += _ion_pill_svg(right_cx, cy, anion_symbol, anion_count, anion_charge)
    svg += "</svg>"
    return svg


# ---------------------------------------------------------------------------
# RDKit 골격 구조식 렌더링 (공유 결합 분자 전용 — 이온성 화합물은 위 알약 카드로 처리)
# ---------------------------------------------------------------------------
FORMULA_SMILES = {normalize_key(f): smiles for f, _, smiles in REAGENTS}


def render_structure_svg(formula: str, width: int = 220, height: int = 150):
    key = normalize_key(formula)

    ion_pair = IONIC_PAIRS.get(key)
    if ion_pair:
        return build_ionic_svg(*ion_pair)

    smiles = FORMULA_SMILES.get(key)
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.drawOptions().clearBackground = False
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
    except Exception:
        return None
    return re.sub(r"<\?xml[^>]*\?>\s*", "", svg).strip()


RESULT_TEMPLATE = (Path(__file__).parent / "result_card.html").read_text(encoding="utf-8")


def build_result_html(result: dict) -> str:
    structure_block = ""
    if result["svg"]:
        structure_block = (
            '<div class="structure-wrap"><div class="structure-svg-holder">'
            + result["svg"]
            + "</div></div>"
        )
    html = RESULT_TEMPLATE
    html = html.replace("__BASIS__", result["basis"])
    html = html.replace("__STRUCTURE_BLOCK__", structure_block)
    html = html.replace("__MOL_FINAL__", repr(float(result["mol"])))
    html = html.replace("__MASS_FINAL__", repr(float(result["mass"])))
    html = html.replace("__PILE_TOP__", result["pile_top"])
    html = html.replace("__PILE_BOTTOM__", result["pile_bottom"])
    return html


# ---------------------------------------------------------------------------
# 화면 구성
# ---------------------------------------------------------------------------
GLOBAL_CSS = """
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Baloo+2:wght@700;800&display=swap" rel="stylesheet">
<style>
  :root {
    --blue-900: #223a63;
    --blue-700: #3a5fad;
    --blue-500: #5478c9;
    --blue-100: #e7e9ff;
    --blue-050: #f1f2ff;
    --ink: #1c2942;
    --muted: #6b7a99;
    --border: #d9defa;
    --error: #e0483e;
    --white: #ffffff;
  }

  .stApp {
    position: relative;
    overflow-x: hidden;
    background: linear-gradient(180deg, #eef0ff 0%, #cdd9ff 100%);
    font-family: "Pretendard", "Malgun Gothic", "Apple SD Gothic Neo", -apple-system,
      BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  }

  #MainMenu, footer, header[data-testid="stHeader"] { visibility: hidden; }

  .block-container {
    position: relative;
    z-index: 1;
    max-width: 560px;
    margin: 0 auto;
    padding-top: 2.2rem;
    padding-bottom: 3rem;
  }

  .bg-deco {
    position: absolute;
    pointer-events: none;
    color: var(--blue-700);
    opacity: 0.13;
    z-index: 0;
  }

  .bg-deco svg { display: block; width: 100%; height: auto; }
  .bg-deco-1 { top: -40px; right: -70px; width: 300px; transform: rotate(18deg); }
  .bg-deco-2 { bottom: -30px; left: -90px; width: 420px; transform: rotate(-14deg); }

  @keyframes flaskFloat {
    0%, 100% { transform: rotate(14deg) translateY(0); }
    50% { transform: rotate(14deg) translateY(-9px); }
  }

  .title-row {
    display: flex;
    align-items: center;
    gap: 0;
    margin-bottom: 24px;
  }

  .flask-icon {
    width: 96px;
    height: 96px;
    flex-shrink: 0;
    margin-left: -20px;
    margin-top: -6px;
    overflow: visible;
    animation: flaskFloat 3.2s ease-in-out infinite;
  }

  .title-row h1 {
    font-family: "Baloo 2", "Jua", "Malgun Gothic", "Apple SD Gothic Neo", sans-serif;
    font-size: 38px;
    font-weight: 800;
    color: var(--blue-700);
    letter-spacing: -0.01em;
    margin: 0;
    line-height: 1.15;
  }

  [data-testid="stWidgetLabel"] p {
    font-size: 13px;
    font-weight: 700;
    color: var(--blue-900);
  }

  div[data-testid="stForm"] {
    background: var(--white);
    border: 1px solid var(--border);
    border-radius: 20px;
    padding: 24px 20px 28px;
    box-shadow: 0 8px 30px rgba(20, 80, 201, 0.08);
  }

  div[data-testid="stTextInput"] input,
  div[data-testid="stNumberInput"] input {
    border: 1.5px solid var(--border) !important;
    border-radius: 12px !important;
    background: var(--blue-050) !important;
    padding: 10px 14px !important;
    font-size: 16px !important;
    color: var(--ink) !important;
  }

  div[data-testid="stTextInput"] input:focus,
  div[data-testid="stNumberInput"] input:focus {
    border-color: var(--blue-500) !important;
    box-shadow: 0 0 0 3px rgba(47, 111, 237, 0.15) !important;
  }

  div[data-testid="stSelectbox"] div[data-baseweb="select"] > div {
    border: 1.5px solid var(--border) !important;
    border-radius: 12px !important;
    background: var(--blue-050) !important;
  }

  div[data-testid="stFormSubmitButton"] {
    display: flex;
    justify-content: center;
    margin-top: 8px;
  }

  div[data-testid="stFormSubmitButton"] button {
    background: linear-gradient(135deg, var(--blue-700), var(--blue-500));
    color: var(--white);
    border: none;
    border-radius: 999px;
    padding: 10px 40px;
    font-size: 15px;
    font-weight: 800;
    letter-spacing: 0.02em;
    box-shadow: 0 6px 16px rgba(58, 95, 173, 0.35);
  }

  div[data-testid="stFormSubmitButton"] button:hover {
    transform: translateY(-1px);
  }

  .calc-error {
    color: var(--error);
    font-size: 13px;
    text-align: center;
    margin: 10px 0 0;
  }

  .app-footer {
    text-align: center;
    color: var(--muted);
    font-size: 12px;
    margin-top: 18px;
  }

  @media (max-width: 480px) {
    .bg-deco-1 { width: 210px; top: -30px; right: -55px; }
    .bg-deco-2 { width: 280px; bottom: -20px; left: -60px; }
    .title-row { position: relative; justify-content: center; }
    .title-row h1 { font-size: 28px; text-align: center; }
    .flask-icon {
      position: absolute;
      top: -2px;
      left: 50%;
      margin-left: 62px;
      margin-top: 0;
      width: 74px;
      height: 74px;
    }
    div[data-testid="stForm"] { padding: 20px 16px 24px; }
  }
</style>
"""

BACKGROUND_HTML = """
<div class="bg-deco bg-deco-1" aria-hidden="true">
  <svg viewBox="0 0 200 200" fill="none" xmlns="http://www.w3.org/2000/svg">
    <polygon points="100,20 169.3,60 169.3,140 100,180 30.7,140 30.7,60"
      stroke="currentColor" stroke-width="4" stroke-linejoin="round"/>
    <line x1="113.86" y1="40" x2="145.05" y2="58" stroke="currentColor" stroke-width="4" stroke-linecap="round"/>
    <line x1="145.05" y1="142" x2="113.86" y2="160" stroke="currentColor" stroke-width="4" stroke-linecap="round"/>
    <line x1="41.1" y1="118" x2="41.1" y2="82" stroke="currentColor" stroke-width="4" stroke-linecap="round"/>
  </svg>
</div>
<div class="bg-deco bg-deco-2" aria-hidden="true">
  <svg viewBox="0 0 340 200" fill="none" xmlns="http://www.w3.org/2000/svg">
    <polygon points="100,20 30.7,60 30.7,140 100,180 169.3,140 169.3,60"
      stroke="currentColor" stroke-width="5" stroke-linejoin="round"/>
    <polygon points="169.3,60 238.6,20 307.9,60 307.9,140 238.6,180 169.3,140"
      stroke="currentColor" stroke-width="5" stroke-linejoin="round"/>
    <line x1="113.86" y1="40" x2="145.05" y2="58" stroke="currentColor" stroke-width="5" stroke-linecap="round"/>
    <line x1="145.05" y1="142" x2="113.86" y2="160" stroke="currentColor" stroke-width="5" stroke-linecap="round"/>
    <line x1="41.1" y1="118" x2="41.1" y2="82" stroke="currentColor" stroke-width="5" stroke-linecap="round"/>
    <line x1="224.74" y1="40" x2="193.55" y2="58" stroke="currentColor" stroke-width="5" stroke-linecap="round"/>
    <line x1="193.55" y1="142" x2="224.74" y2="160" stroke="currentColor" stroke-width="5" stroke-linecap="round"/>
    <line x1="297.5" y1="118" x2="297.5" y2="82" stroke="currentColor" stroke-width="5" stroke-linecap="round"/>
  </svg>
</div>
"""

HEADER_HTML = """
<div class="title-row">
  <svg class="flask-icon" viewBox="0 0 64 64" xmlns="http://www.w3.org/2000/svg">
    <defs>
      <linearGradient id="glassGrad" x1="20" y1="8" x2="52" y2="58" gradientUnits="userSpaceOnUse">
        <stop offset="0%" stop-color="#ffffff" stop-opacity="0.85"/>
        <stop offset="55%" stop-color="#d9def0" stop-opacity="0.55"/>
        <stop offset="100%" stop-color="#a8b6dd" stop-opacity="0.45"/>
      </linearGradient>
      <linearGradient id="liquidGrad" x1="14" y1="34" x2="50" y2="58" gradientUnits="userSpaceOnUse">
        <stop offset="0%" stop-color="#7fa0d9"/>
        <stop offset="55%" stop-color="#5478c9"/>
        <stop offset="100%" stop-color="#2c4270"/>
      </linearGradient>
      <linearGradient id="neckGrad" x1="26" y1="6" x2="38" y2="6" gradientUnits="userSpaceOnUse">
        <stop offset="0%" stop-color="#a3aed4"/>
        <stop offset="50%" stop-color="#f2f7ff"/>
        <stop offset="100%" stop-color="#a3aed4"/>
      </linearGradient>
      <filter id="flaskShadow" x="-40%" y="-20%" width="180%" height="160%">
        <feDropShadow dx="0" dy="3" stdDeviation="2.4" flood-color="#3a5fad" flood-opacity="0.35"/>
      </filter>
    </defs>
    <g filter="url(#flaskShadow)">
      <path d="M28 7V24.5L14.5 48.5C12.7 51.7 15 55.7 18.7 55.7H45.3C49 55.7 51.3 51.7 49.5 48.5L36 24.5V7Z"
        fill="url(#glassGrad)" stroke="#3a5fad" stroke-width="2.3" stroke-linejoin="round"/>
      <path d="M16.7 44.8C15.6 46.8 14.6 48.9 15.6 51.2C16.6 53.4 18.9 53.9 21.3 53.9H43.4C45.6 53.9 47.7 53.3 48.5 51.3C49.4 49.1 48.4 46.8 47.2 44.8C43.6 46.6 38.3 47.7 32 47.7C25.7 47.7 20.4 46.6 16.7 44.8Z"
        fill="url(#liquidGrad)"/>
      <ellipse cx="24" cy="49.5" rx="3.4" ry="1.6" fill="#ffffff" opacity="0.35"/>
      <circle cx="34" cy="51" r="2" fill="#bfe0ff" opacity="0.9"/>
      <circle cx="40" cy="44.5" r="1.4" fill="#eaf3ff" opacity="0.9"/>
      <circle cx="29" cy="43.5" r="1.1" fill="#eaf3ff" opacity="0.8"/>
      <path d="M29.2 9V24.5" stroke="#ffffff" stroke-width="1.6" stroke-linecap="round" opacity="0.7"/>
      <rect x="26" y="5.2" width="12" height="3.4" rx="1.6" fill="url(#neckGrad)" stroke="#3a5fad" stroke-width="1.6"/>
      <circle cx="30" cy="4" r="1.6" fill="#8fc0ff" opacity="0">
        <animate attributeName="cy" values="4;-15" dur="2.4s" begin="0s" repeatCount="indefinite"/>
        <animate attributeName="cx" values="30;27;30" dur="2.4s" begin="0s" repeatCount="indefinite"/>
        <animate attributeName="opacity" values="0;0.9;0" keyTimes="0;0.3;1" dur="2.4s" begin="0s" repeatCount="indefinite"/>
        <animate attributeName="r" values="1;2;1.2" dur="2.4s" begin="0s" repeatCount="indefinite"/>
      </circle>
      <circle cx="33" cy="4" r="1.2" fill="#bfe0ff" opacity="0">
        <animate attributeName="cy" values="4;-13" dur="2s" begin="0.6s" repeatCount="indefinite"/>
        <animate attributeName="cx" values="33;36;33" dur="2s" begin="0.6s" repeatCount="indefinite"/>
        <animate attributeName="opacity" values="0;0.9;0" keyTimes="0;0.3;1" dur="2s" begin="0.6s" repeatCount="indefinite"/>
        <animate attributeName="r" values="0.8;1.6;1" dur="2s" begin="0.6s" repeatCount="indefinite"/>
      </circle>
      <circle cx="31.5" cy="4" r="1" fill="#eaf3ff" opacity="0">
        <animate attributeName="cy" values="4;-17" dur="2.8s" begin="1.2s" repeatCount="indefinite"/>
        <animate attributeName="cx" values="31.5;29;31.5" dur="2.8s" begin="1.2s" repeatCount="indefinite"/>
        <animate attributeName="opacity" values="0;0.85;0" keyTimes="0;0.3;1" dur="2.8s" begin="1.2s" repeatCount="indefinite"/>
        <animate attributeName="r" values="0.6;1.4;0.8" dur="2.8s" begin="1.2s" repeatCount="indefinite"/>
      </circle>
    </g>
  </svg>
  <h1>Experimental Helper</h1>
</div>
"""

st.html(GLOBAL_CSS)
st.html(BACKGROUND_HTML)
st.html(HEADER_HTML)

quick_options = ["직접 입력"] + [f"{f}  ·  {n}" for f, n in COMMON_REAGENTS]


def _apply_quick_pick():
    choice = st.session_state.get("quick_pick")
    if choice and choice != "직접 입력":
        st.session_state["formula_text"] = choice.split("  ·  ", 1)[0]


st.selectbox("자주 쓰는 시약", quick_options, key="quick_pick", on_change=_apply_quick_pick)

with st.form("calc_form"):
    formula = st.text_input("화학식", key="formula_text", placeholder="예: NaCl, CuSO4·5H2O")
    col1, col2 = st.columns(2)
    with col1:
        molarity = st.number_input("몰농도 (M)", min_value=0.0, value=0.0, step=0.01, format="%.4f", key="molarity_input")
    with col2:
        volume = st.number_input("최종 부피 (mL)", min_value=0.0, value=0.0, step=1.0, format="%.2f", key="volume_input")
    submitted = st.form_submit_button("계산")

if submitted:
    error = None
    mass = None

    if not formula.strip():
        error = "화학식을 입력해주세요"
    else:
        try:
            mass = molar_mass(formula)
        except Exception:
            error = "화학식을 확인해주세요 (예: CuSO4·5H2O)"

    if not error and molarity <= 0:
        error = "몰농도를 입력해주세요"
    if not error and volume <= 0:
        error = "최종 부피를 입력해주세요"

    if error:
        st.session_state["calc_error"] = error
        st.session_state["calc_result"] = None
    else:
        volume_l = volume / 1000
        final_mol = molarity * volume_l
        final_mass = final_mol * mass
        pile_top, pile_bottom = get_pile_colors(formula)
        st.session_state["calc_error"] = None
        st.session_state["calc_result"] = {
            "basis": f"{format_number(volume)} mL 기준",
            "mol": final_mol,
            "mass": final_mass,
            "svg": render_structure_svg(formula),
            "pile_top": pile_top,
            "pile_bottom": pile_bottom,
        }

if st.session_state.get("calc_error"):
    st.html(f'<p class="calc-error">{st.session_state["calc_error"]}</p>')

result = st.session_state.get("calc_result")
if result:
    height = 470 if result["svg"] else 340
    components.html(build_result_html(result), height=height, scrolling=False)

st.html('<div class="app-footer">몰질량(M) × 몰농도(mol/L) × 부피(L) = 무게(g)</div>')
