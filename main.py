import streamlit as st
from data.multiapp import MultiApp
from data.ot2_cloning import app_v2 as ot2_cloning
st.set_page_config(layout="wide")
app = MultiApp()

# Add all your application here
app.add_app("OT-2 Cloning", ot2_cloning.main)

# The main app
if __name__ == "__main__":
    app.run()
