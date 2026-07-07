from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(page_title="Experimental Helper", page_icon="🧪", layout="wide")

st.markdown(
    """
    <style>
        #MainMenu, header, footer { visibility: hidden; }
        .block-container { padding: 0 !important; max-width: 100% !important; }
        iframe { border: none; display: block; }
    </style>
    """,
    unsafe_allow_html=True,
)

html_content = (Path(__file__).parent / "index.html").read_text(encoding="utf-8")
components.html(html_content, height=1000, scrolling=True)
