import streamlit as st
import os

st.set_page_config(
    page_title="TAILVAR: Terminal codon Analysis and Improved prediction of Lengthened VARiants",
    page_icon="🧬",
    layout="centered"
)

# Title
st.title("🧬 TAILVAR")
st.markdown("""
Welcome to the TAILVAR website to analyze pathogenicity predictions for **stop-loss variants**

Upload the variant information of your interest and receive the predicted TAILVAR score.
""")

# Sidebar
with st.sidebar:
    st.header("About TAILVAR")
    st.markdown("""
    - Developed for stop-loss variant analysis
    - Supports `.vcf` or `.txt` formats
    - Results include TAILVAR scores and classifications
    """)
    st.markdown("[🔗 GitHub Repository](https://github.com/dr-toon/TAILVAR)")
    st.markdown("Contact: YOONJH@yuhs.ac")

# Upload file
st.subheader("📁 Upload Your File")

uploaded_file = st.file_uploader(
    "Upload a VCF or TXT file (max 2MB)", type=["vcf", "txt"]
)

if uploaded_file:
    upload_path = os.path.join("uploads", uploaded_file.name)
    os.makedirs("uploads", exist_ok=True)

    with open(upload_path, "wb") as f:
        f.write(uploaded_file.read())

    st.success(f"✅ File '{uploaded_file.name}' uploaded successfully!")

    # Analysis (not functional yet)
    if st.button("🚀 Get TAILVAR score"):
        st.info("📡 Running analysis... (to be implemented)")
        # 이후 분석 연동: subprocess로 bash script 호출

else:
    st.warning("⬆️ Please upload a file to begin.")

# Footer
st.markdown("---")
st.caption("Developed by Jihoon Yoon, M.D., Ph.D. | Powered by Streamlit")
