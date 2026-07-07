import re
import streamlit as st

st.set_page_config(page_title="Experimental Helper", page_icon="🧪", layout="centered")

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

COMMON_REAGENTS = [
    ("NaCl", "염화나트륨"),
    ("NaOH", "수산화나트륨"),
    ("KOH", "수산화칼륨"),
    ("HCl", "염산"),
    ("H2SO4", "황산"),
    ("NaHCO3", "탄산수소나트륨"),
    ("Na2CO3", "탄산나트륨"),
    ("K2CO3", "탄산칼륨"),
    ("CaCl2", "염화칼슘"),
    ("KCl", "염화칼륨"),
    ("NH4Cl", "염화암모늄"),
    ("NaNO3", "질산나트륨"),
    ("KNO3", "질산칼륨"),
    ("AgNO3", "질산은"),
    ("FeCl3", "염화철(III)"),
    ("CuSO4·5H2O", "황산구리(II) 오수화물"),
    ("MgSO4·7H2O", "황산마그네슘 칠수화물"),
    ("Na2SO4", "황산나트륨"),
    ("Ca(OH)2", "수산화칼슘"),
    ("CH3COOH", "아세트산"),
    ("C2H5OH", "에탄올"),
    ("C6H12O6", "포도당"),
]

STRUCTURE_MAP = {
    "H2O": "💧",
    "HCL": "🧪",
    "NH3": "☁️",
    "CO2": "🌫️",
    "NACL": "🧂",
    "KCL": "🧂",
    "NAOH": "🧴",
    "KOH": "🧴",
    "CH3COOH": "🧴",
    "C2H5OH": "🧴",
}


def parse_formula(formula: str) -> dict:
    parts = [part.strip() for part in re.split(r"[·.\*]", formula) if part.strip()]
    if not parts:
        raise ValueError("empty")

    total = {}
    for part in parts:
        match = re.match(r"^(\d+)?(.+)$", part)
        if not match:
            raise ValueError("invalid")
        multiplier = int(match.group(1) or 1)
        body = match.group(2)
        counts = parse_group(body)
        for element, count in counts.items():
            total[element] = total.get(element, 0) + count * multiplier
    return total


def parse_group(formula: str) -> dict:
    pos = 0

    def read_number() -> int:
        nonlocal pos
        start = pos
        while pos < len(formula) and formula[pos].isdigit():
            pos += 1
        return int(formula[start:pos]) if pos > start else 1

    def parse_expr() -> dict:
        nonlocal pos
        counts = {}
        while pos < len(formula):
            ch = formula[pos]
            if ch in "([":
                close = ")" if ch == "(" else "]"
                pos += 1
                inner = parse_expr()
                if pos >= len(formula) or formula[pos] != close:
                    raise ValueError("mismatched brackets")
                pos += 1
                multiplier = read_number()
                for element, count in inner.items():
                    counts[element] = counts.get(element, 0) + count * multiplier
            elif ch in ")]":
                break
            elif ch.isupper():
                element = ch
                pos += 1
                if pos < len(formula) and formula[pos].islower():
                    element += formula[pos]
                    pos += 1
                if element not in ATOMIC_WEIGHTS:
                    raise ValueError(f"unknown element: {element}")
                multiplier = read_number()
                counts[element] = counts.get(element, 0) + multiplier
            else:
                raise ValueError(f"unexpected character: {ch}")
        return counts

    result = parse_expr()
    if pos != len(formula):
        raise ValueError("unexpected character")
    return result


def molar_mass(formula: str) -> float:
    counts = parse_formula(formula)
    return sum(ATOMIC_WEIGHTS[element] * count for element, count in counts.items())


def format_number(value: float, digits: int = 4) -> str:
    if not value or not float(value):
        return "0"
    return f"{value:,.{digits}f}"


def render_structure(formula: str) -> str:
    normalized = formula.replace(" ", "").upper()
    return STRUCTURE_MAP.get(normalized, "⚗️")


st.markdown(
    """
    <style>
        .stApp { background: linear-gradient(135deg, #eef2ff 0%, #dfe8ff 100%); }
        .block-container { padding-top: 2rem; padding-bottom: 2rem; }
        div[data-testid="stMetric"] { background: white; border: 1px solid #d9e2ff; border-radius: 14px; padding: 12px; }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("Experimental Helper")
st.caption("몰질량 기반으로 필요한 시약의 몰수와 질량을 계산합니다.")

with st.form("calc_form"):
    formula = st.text_input(
        "화학식",
        value="",
        placeholder="예: NaCl, CuSO4·5H2O",
        help="분자량이 자동으로 계산됩니다.",
    )
    col1, col2 = st.columns(2)
    with col1:
        molarity = st.number_input("몰농도 (M)", min_value=0.0, step=0.01, format="%.6f")
    with col2:
        volume_ml = st.number_input("최종 부피 (mL)", min_value=0.0, step=0.1, format="%.6f")

    selected_formula = st.selectbox("자주 쓰는 시약", ["직접 입력"] + [f"{formula_name} ({name})" for formula_name, name in COMMON_REAGENTS], index=0)
    submitted = st.form_submit_button("계산")

    if selected_formula != "직접 입력":
        selected_formula_name = selected_formula.split(" (")[0]
        formula = selected_formula_name

if submitted:
    if not formula.strip():
        st.error("화학식을 입력해주세요.")
    elif molarity <= 0:
        st.error("몰농도를 입력해주세요.")
    elif volume_ml <= 0:
        st.error("최종 부피를 입력해주세요.")
    else:
        try:
            mass_per_mole = molar_mass(formula)
            volume_l = volume_ml / 1000.0
            moles = molarity * volume_l
            mass_g = moles * mass_per_mole
        except Exception:
            st.error("화학식을 확인해주세요. 예: CuSO4·5H2O")
        else:
            st.success(f"몰질량: {format_number(mass_per_mole)} g/mol")
            st.write(f"{volume_ml} mL 기준")
            st.metric("넣어야 하는 시약의 몰수", f"{moles:.6f} mol")
            st.metric("넣어야 하는 시약의 질량", f"{mass_g:.6f} g")
            st.caption(f"계산식: 몰질량 × 몰농도 × 부피(L) = 무게(g)")
            st.markdown(f"**구조/형태 예시:** {render_structure(formula)}")

elif formula.strip():
    try:
        mass_per_mole = molar_mass(formula)
    except Exception:
        st.info("화학식을 확인해주세요.")
    else:
        st.caption(f"몰질량: {format_number(mass_per_mole)} g/mol")
