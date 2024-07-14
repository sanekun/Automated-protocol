import streamlit as st
import pandas as pd
import re
from io import StringIO
import json
from datetime import datetime
from data.ot2.protocols.cloning.check_protocol import *

# Statics
def empty_plate_df():
    # Empty 96-well plate DataFrame
    rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
    columns = [str(i) for i in range(1, 13)]
    plate_df = pd.DataFrame(index=rows, columns=columns)
    return plate_df


def plate_initialization(i, plate_type: str, table_type="empty"):
    # Make Default Plate Map
    # Make Default session state value
    initial_tables = {
        "empty": empty_plate_df(),
        "PCR": pd.DataFrame(
            {
                "Name": [None for _ in range(1)],
                "DNA1": [None for _ in range(1)],
                "DNA2": [None for _ in range(1)],
                "DNA3": [None for _ in range(1)],
                "Enzyme1": ["PCRMix" for _ in range(1)],
                "DW": ["DW" for _ in range(1)],
            }
        ),
        "Assembly": pd.DataFrame(
            {
                "Name": [None for _ in range(1)],
                "DNA1": [None for _ in range(1)],
                "DNA2": [None for _ in range(1)],
                "Enzyme1": ["AssemblyMix" for _ in range(1)],
                "DW": ["DW" for _ in range(1)],
            }
        ),
    }
    if table_type == "empty":
        if f"{plate_type}_plate_{i}_name" not in st.session_state:
            st.session_state[
                f"{plate_type}_plate_{i}_name"
            ] = f"{plate_type}_plate_{i}"
        if f"{plate_type}_plate_{i}_toggle" not in st.session_state:
            st.session_state[f"{plate_type}_plate_{i}_toggle"] = False

    if f"{plate_type}_plate_{i}_df" not in st.session_state:
        st.session_state[f"{plate_type}_plate_{i}_df"] = initial_tables[table_type]

def plate_transformation(df, data_form):
    # Change data form to [long or wide]
    assert data_form in ["wide", "long"], "Plate_transformation: data_form Error"
    if data_form == "long":
        new_df = (
            df.reset_index()
            .melt(id_vars=["index"], value_name="Value")
            .rename(columns={"index": "Row", "variable": "Column"})
        )
        new_df["well"] = new_df["Row"] + new_df["Column"].astype(str)
        new_df = new_df.set_index("well")["Value"].to_frame()

        return new_df

    elif data_form == "wide":
        new_df = df.reset_index()
        new_df[["Row", "Column"]] = new_df["well"].str.extract(r"([A-Z]+)(\d+)")
        new_df["Column"] = new_df["Column"].astype(int)
        new_df = new_df[["Row", "Column", "Value"]]

        new_df = new_df.pivot(index="Row", columns="Column", values="Value").sort_index(
            axis=1
        )
        new_df.columns.name = None
        new_df.index.name = None

        return new_df


def editor_update(editor_key, df, df_form):
    # For convert long form to wide form

    # Cuz of the st.data_editor's behavior.
    # It store the edited data in the dictionary
    # And if script re-run by another behavior If update the dataframe and send to output.
    # If I make the output will be the original dataframe, It updated every time when i change table and re-read again for make data_editor.
    # If makes some errors to editing tables.

    # So if we change toogle value, we need to update directly the dataframe and transform.
    if df_form == "wide":
        editor = st.session_state[editor_key]["edited_rows"]
        for rows, updates in editor.items():
            for cols, value in updates.items():
                df.loc[df.index[rows], cols] = value
        editor = st.session_state[editor_key]["added_rows"]
        for i in editor:
            df.loc[len(df)] = i

    if df_form == "long":
        editor = st.session_state[editor_key]["edited_rows"]
        for rows, updates in editor.items():
            row, col = re.match(
                "([A-Z]+)(\d+)", plate_transformation(df, "long").index[rows]
            ).groups()
            for _, value in updates.items():
                df.loc[row, col] = value


def add_column_button(df, column_type):
    # DNA, Enzyme 등 마지막 열 옆에 새로운 열 추가

    # Check current df's column name (like DNA1, DNA2, DNA3) add final number
    column_names = df.columns
    # grep column_type+[0-9]+ string in column_names
    column_numbers = [
        int(re.findall(rf"{column_type}(\d+)", column_name)[0])
        for column_name in column_names
        if re.findall(rf"{column_type}(\d+)", column_name)
    ]
    df[f"{column_type}{max(column_numbers)+1}"] = [None for _ in range(len(df))]
    # column order is 'well', 'name', ['DNA\d+', 'Enzyme[+]', 'DW'']
    df = df.reindex(
        ["Name"]
        + sorted(
            [
                column_name
                for column_name in df.columns
                if re.findall(rf"DNA(\d+)", column_name)
            ]
        )
        + sorted(
            [
                column_name
                for column_name in df.columns
                if re.findall(rf"Enzyme(\d+)", column_name)
            ]
        )
        + ["DW"],
        axis=1,
    )
    return df


def check_overlap_name(df1, df2):
    # Reaction, TF table에서 등록한 Name만 활용하기 위하여 Overlap 이름 확인
    # Check overlap in each dataframe first If overlap name exist, return False
    # Check overlap between df1 and df2
    # If overlap name exist, return False
    df1_names = df1.dropna().unique()
    df2_names = df2.dropna().unique()
    if len(df1_names) != len(set(df1_names)):
        return False
    if len(df2_names) != len(set(df2_names)):
        return False
    if len(set(df1_names) & set(df2_names)) != 0:
        return False
    return list(df1_names) + list(df2_names)

@st.cache_data
def read_README(path):
    README = []
    with open(path, "r") as f:
        txt=f.readlines()

    txt.index("# Developer only\n")
    README.append("".join(txt[:txt.index("# Developer only\n")]))
    README.append("".join(txt[txt.index("# Developer only\n"):]))

    return README

def Use_example():
    # Example Data set
    # DNA table
    st.session_state["source_plate_0_df"].loc["A", "1"] = "Template1"
    st.session_state["source_plate_0_df"].loc["B", "1"] = "Primer1"
    st.session_state["source_plate_0_df"].loc["C", "1"] = "Primer2"

    st.session_state["source_plate_0_df"].loc["A", "2"] = "Template2"
    st.session_state["source_plate_0_df"].loc["B", "2"] = "Primer3"
    st.session_state["source_plate_0_df"].loc["C", "2"] = "Primer4"

    # PCR & Assembly table
    st.session_state["PCR_plate_0_df"].loc[0, "Name"] = "Vector"
    st.session_state["PCR_plate_0_df"].loc[0, "DNA1"] = "Template1"
    st.session_state["PCR_plate_0_df"].loc[0, "DNA2"] = "Primer1"
    st.session_state["PCR_plate_0_df"].loc[0, "DNA3"] = "Primer2"
    st.session_state["PCR_plate_0_df"].loc[0, "Enzyme1"] = "PCRMix"
    st.session_state["PCR_plate_0_df"].loc[0, "DW"] = "DW"

    st.session_state["PCR_plate_0_df"].loc[1, "Name"] = "Insert"
    st.session_state["PCR_plate_0_df"].loc[1, "DNA1"] = "Template2"
    st.session_state["PCR_plate_0_df"].loc[1, "DNA2"] = "Primer3"
    st.session_state["PCR_plate_0_df"].loc[1, "DNA3"] = "Primer4"
    st.session_state["PCR_plate_0_df"].loc[1, "Enzyme1"] = "PCRMix"
    st.session_state["PCR_plate_0_df"].loc[1, "DW"] = "DW"

    st.session_state["Assembly_plate_0_df"].loc[0, "Name"] = "Assembled_vector"
    st.session_state["Assembly_plate_0_df"].loc[0, "DNA1"] = "Vector"
    st.session_state["Assembly_plate_0_df"].loc[0, "DNA2"] = "Insert"

    # Reaction table
    st.session_state["Reaction_plate_0_df"].loc["A", "1"] = "Vector"
    st.session_state["Reaction_plate_0_df"].loc["B", "1"] = "Insert"
    st.session_state["Reaction_plate_0_df"].loc["A", "3"] = "Assembled_vector"
    st.session_state["TF_plate_0_df"].loc["A", "1"] = "Assembled_vector"
    st.session_state["TF_plate_0_df"].loc["A", "2"] = "Assembled_vector"
    st.session_state["TF_plate_0_df"].loc["A", "3"] = "Assembled_vector"


def load_previous_protocol():
    if st.session_state['load'] is None:
        return None
    
    # Read and Convert to String
    data = StringIO(st.session_state['load'].getvalue().decode("utf-8")).readlines()
    idx = data.index("# [PARAMETERS]\n")
    params = data[idx+1]
    params = params.replace("PARAMETERS = ", "")
    params = eval(params)

    if params['Parameters']['protocol'] != "OT-2 cloning":
        st.error("This is not OT-2 cloning protocol")
        return None

    # Fill parameters
    # st.session_state['messenger'] = params['Parameters']['Messenger']
    st.session_state['stop_between_reactions'] = params['Parameters']['Stop_between_reactions']
    st.session_state['PCR_extension_time'] = params['Parameters']['PCR_extension_time']
    st.session_state['TF_recovery_time'] = params['Parameters']['TF_recovery_time']

    # Source plate
    for target_type in params["Parameters"]["Plate_type"]:
        target_plates = [plate for plate in params["Plates"].values() if plate["type"] == target_type]
        st.session_state[f"num_of_{target_type}_plate"] = len(target_plates)
        
        for n, plate in enumerate(target_plates):
            st.session_state[f"{target_type}_plate_{n}_name"] = plate["name"]
            st.session_state[f"{target_type}_plate_{n}_labware"] = plate["labware"]
            st.session_state[f"{target_type}_plate_{n}_toggle"] = False
            
            # Read plate data and Update to existed dataframe
            st.session_state[f"{target_type}_plate_{n}_df"] = empty_plate_df()
            tmp = pd.DataFrame.from_dict(plate["data"], orient="index", columns=['Value'])
            tmp.index.name = 'well'
            tmp2 = plate_transformation(tmp, 'wide')
            tmp2.columns = tmp2.columns.map(str)
            st.session_state[f"{target_type}_plate_{n}_df"].update(tmp2)
    
    # Reaction Table & Volumes
    for reaction in params["Reactions"].values():
        st.session_state[f"{reaction['type']}_plate_{0}_df"] = pd.DataFrame.from_dict(reaction['data'])
        # st.session_state[f"{reaction['type']}_volume"] = pd.DataFrame.from_dict(params["Reaction_volume"][f"{reaction}"], orient='index').T
        # Volume 입력 구조를 변경해야함.

def main():
    plate_types = ["source", "Reaction", "TF"]
    for i in plate_types:
        if f"num_of_{i}_plate" not in st.session_state:
            st.session_state[f"num_of_{i}_plate"] = 1
    
    st.session_state["labwares"] = ["biorad_96_wellplate_200ul_pcr", "nest_96_wellplate_200ul_flat"]
    export_JSON = False
    README = read_README('data/ot2/protocols/cloning/README.md')

    # Main
    with st.expander("Manual", expanded=True):
        manual_tabs = st.tabs(["README", "Developer only"])
        with manual_tabs[0]:
            with st.container(height=450):
                st.markdown(
                    README[0]
                )

        with manual_tabs[1]:
            with st.container(height=450):
                st.markdown(
                    README[1]
                )

    st.button("Use Example", on_click=Use_example)
    st.file_uploader("## Load previous result", type=['py'],key="load", disabled=False, help="Load previous result to modify", on_change=load_previous_protocol)
    st.markdown("---")

    ## DNA Plate
    st.markdown("# Main")
    st.markdown("## Source Plate")
    with st.expander("Source plate", expanded=True):
        plate_type = "source"
        number_of_plate = st.number_input(
            f"Number of {plate_type}_plate",
            min_value=1,
            step=1,
            value=st.session_state[f"num_of_{plate_type}_plate"],
            key=f"num_of_{plate_type}_plate",
        )
        plates = st.tabs([f"{plate_type}_Plate_{i}" for i in range(number_of_plate)])

        for i in range(len(plates)):
            plate_initialization(i, plate_type=plate_type, table_type="empty")
            with plates[i]:
                st.text_input(
                    "Name",
                    placeholder=f"{plate_type}_plate_{i}_name",
                    key=f"{plate_type}_plate_{i}_name",
                )
                st.selectbox(
                    "Labware",
                    options=st.session_state["labwares"],
                    key=f"{plate_type}_plate_{i}_labware",
                )

                st.toggle(
                    "Long Form",
                    key=f"{plate_type}_plate_{i}_toggle",
                    on_change=editor_update, # Transformation을 여기에서 할걸
                    kwargs={
                        "editor_key": f"{plate_type}_plate_{i}_long"
                        if st.session_state[f"{plate_type}_plate_{i}_toggle"]
                        else f"{plate_type}_plate_{i}_wide",
                        "df": st.session_state[f"{plate_type}_plate_{i}_df"],
                        "df_form": "long"
                        if st.session_state[f"{plate_type}_plate_{i}_toggle"]
                        else "wide",
                    },
                )
                if not st.session_state[f"{plate_type}_plate_{i}_toggle"]:
                    st.data_editor(
                        data=st.session_state[f"{plate_type}_plate_{i}_df"],
                        use_container_width=True,
                        key=f"{plate_type}_plate_{i}_wide",
                    )
                else:
                    st.data_editor(
                        data=plate_transformation(
                            st.session_state[f"{plate_type}_plate_{i}_df"], data_form="long"
                        ),
                        use_container_width=True,
                        key=f"{plate_type}_plate_{i}_long",
                    )
                9
    st.markdown("## Reaction")
    with st.expander("Reaction - PCR", expanded=True):
        i, reaction_type = 0, "PCR"
        plate_initialization(i, plate_type=reaction_type, table_type="PCR")
        col1, col2 = st.columns([10, 2])
        with col1:
            st.data_editor(
                data=st.session_state[f"{reaction_type}_plate_{i}_df"],
                num_rows="dynamic",
                use_container_width=True,
                key=f"{reaction_type}_plate_{i}_wide",
                hide_index=True,
            )
        with col2:
            with st.container(border=True):
                column_type = st.selectbox(
                    "Add Column",
                    ["DNA", "Enzyme"],
                    label_visibility="collapsed",
                    key=f"{reaction_type}_plate_{i}_selecttype",
                )
                if st.button(
                    "Add Column",
                    key=f"{reaction_type}_plate_{i}_addcolumn",
                    on_click=editor_update,
                    kwargs={
                        "editor_key": f"{reaction_type}_plate_{i}_wide",
                        "df": st.session_state[f"{reaction_type}_plate_{i}_df"],
                        "df_form": "wide",
                    },
                ):
                    st.session_state[f"{reaction_type}_plate_{i}_df"] = add_column_button(
                        df=st.session_state[f"{reaction_type}_plate_{i}_df"],
                        column_type=column_type,
                    )
                    st.rerun()

    with st.expander("Reaction - Assembly", expanded=True):
        i, reaction_type = 0, "Assembly"
        plate_initialization(i, plate_type=reaction_type, table_type="Assembly")
        col1, col2 = st.columns([10, 2])
        with col1:
            st.data_editor(
                data=st.session_state[f"{reaction_type}_plate_{i}_df"],
                num_rows="dynamic",
                use_container_width=True,
                key=f"{reaction_type}_plate_{i}_wide",
                hide_index=True,
            )
        with col2:
            with st.container(border=True):
                column_type = st.selectbox(
                    "Add Column",
                    ["DNA", "Enzyme"],
                    label_visibility="collapsed",
                    key=f"{reaction_type}_plate_{i}_selecttype",
                )
                if st.button(
                    "Add Column",
                    key=f"{reaction_type}_plate_{i}_addcolumn",
                    on_click=editor_update,
                    kwargs={
                        "editor_key": f"{reaction_type}_plate_{i}_wide",
                        "df": st.session_state[f"{reaction_type}_plate_{i}_df"],
                        "df_form": "wide",
                    },
                ):
                    st.session_state[f"{reaction_type}_plate_{i}_df"] = add_column_button(
                        df=st.session_state[f"{reaction_type}_plate_{i}_df"],
                        column_type=column_type,
                    )
                    st.rerun()

    with st.expander("Reaction Plate", expanded=True):
        plate_type = "Reaction"
        number_of_plate = st.number_input(
            f"Number of {plate_type}_plate",
            min_value=1,
            step=1,
            value=st.session_state[f"num_of_{plate_type}_plate"],
            key=f"num_of_{plate_type}_plate",
            disabled=True,
            help="This Method can only one reaction Plate",
        )
        plates = st.tabs([f"{plate_type}_Plate_{i}" for i in range(number_of_plate)])

        for i in range(len(plates)):
            plate_initialization(i, plate_type=plate_type, table_type="empty")
            with plates[i]:
                st.text_input("Name", key=f"{plate_type}_plate_{i}_name")
                st.selectbox(
                    "Labware",
                    options=st.session_state["labwares"],
                    key=f"{plate_type}_plate_{i}_labware",
                )

                st.toggle(
                    "Long Form",
                    key=f"{plate_type}_plate_{i}_toggle",
                    on_change=editor_update,
                    kwargs={
                        "editor_key": f"{plate_type}_plate_{i}_long"
                        if st.session_state[f"{plate_type}_plate_{i}_toggle"]
                        else f"{plate_type}_plate_{i}_wide",
                        "df": st.session_state[f"{plate_type}_plate_{i}_df"],
                        "df_form": "long"
                        if st.session_state[f"{plate_type}_plate_{i}_toggle"]
                        else "wide",
                    },
                )
                # Wide form
                if not st.session_state[f"{plate_type}_plate_{i}_toggle"]:
                    st.data_editor(
                        data=st.session_state[f"{plate_type}_plate_{i}_df"],
                        use_container_width=True,
                        key=f"{plate_type}_plate_{i}_wide",
                        disabled=["2", "4", "6", "8", "10", "12"]
                    )
                else:
                # Long form
                    st.data_editor(
                        data=plate_transformation(
                            st.session_state[f"{plate_type}_plate_{i}_df"], data_form="long"
                        ),
                        use_container_width=True,
                        key=f"{plate_type}_plate_{i}_long",
                    )
            st.success('짝수 Column은 정제 이유로 사용 불가합니다.')

    ## Transformation
    st.markdown("## Transformation")
    with st.expander("Transformation", expanded=True):
        plate_type = "TF"
        number_of_plate = st.number_input(
            f"Number of {plate_type}_plate",
            min_value=1,
            step=1,
            value=st.session_state[f"num_of_{plate_type}_plate"],
            key=f"num_of_{plate_type}_plate",
        )
        plates = st.tabs([f"{plate_type}_Plate_{i}" for i in range(number_of_plate)])

        for i in range(len(plates)):
            plate_initialization(i, plate_type=plate_type, table_type="empty")
            with plates[i]:
                st.text_input("Name", key=f"{plate_type}_plate_{i}_name")
                st.selectbox(
                    "Labware",
                    options=st.session_state["labwares"],
                    key=f"{plate_type}_plate_{i}_labware",
                )

                st.toggle(
                    "Long Form",
                    key=f"{plate_type}_plate_{i}_toggle",
                    on_change=editor_update,
                    kwargs={
                        "editor_key": f"{plate_type}_plate_{i}_long"
                        if st.session_state[f"{plate_type}_plate_{i}_toggle"]
                        else f"{plate_type}_plate_{i}_wide",
                        "df": st.session_state[f"{plate_type}_plate_{i}_df"],
                        "df_form": "long"
                        if st.session_state[f"{plate_type}_plate_{i}_toggle"]
                        else "wide",
                    },
                )
                if not st.session_state[f"{plate_type}_plate_{i}_toggle"]:
                    st.data_editor(
                        data=st.session_state[f"{plate_type}_plate_{i}_df"],
                        use_container_width=True,
                        key=f"{plate_type}_plate_{i}_wide",
                    )
                else:
                    st.data_editor(
                        data=plate_transformation(
                            st.session_state[f"{plate_type}_plate_{i}_df"], data_form="long"
                        ),
                        use_container_width=True,
                        key=f"{plate_type}_plate_{i}_long",
                    )

    st.markdown("## Advanced")
    with st.expander("Advanced", expanded=False):
        st.success("Please adjust here as last step")
        # Reaction-PCR: DNA, Enzyme, DW(up to)
        # Reaction-Assembly: DNA, Enzyme, DW(up to)

        # Make Empty DataFrame with existing column name of PCR and Assembly
        PCR_volume = pd.DataFrame(columns=st.session_state["PCR_plate_0_df"].columns).drop(
            columns=["Name"]
        )
        PCR_volume.loc[0, ["DNA1", "DNA2", "DNA3", "Enzyme1", "DW"]] = [
            "0.5",
            "0.5",
            "0.5",
            "12.5",
            "21",
        ]
        # Thermocycler condition

        Assembly_volume = pd.DataFrame(
            columns=st.session_state["Assembly_plate_0_df"].columns
        ).drop(columns=["Name"])
        
        Assembly_volume.loc[0, ["DNA1", "DNA2", "Enzyme1", "DW"]] = ["1.5", "1.5", "4", "8"]

        st.markdown("### PCR volume")
        update_PCR_volume = st.data_editor(
            data=PCR_volume, use_container_width=True, key="PCR_volume", hide_index=True
        )
        st.markdown("### Assembly volume")
        update_Assembly_volume = st.data_editor(
            data=Assembly_volume,
            use_container_width=True,
            key="Assembly_volume",
            hide_index=True,
        )
        st.checkbox('Stop between Reactions', value=True,
                    key='stop_between_reactions',
                    help='If you want to stop between reactions for take enzymes, check this box')

        advanced_column = st.columns([5,2,5])
        with advanced_column[0]:
            with st.container(border=True):
                st.number_input("PCR extension time (seconds)", min_value=1, step=1, value=15,
                                key='PCR_extension_time')
                st.number_input("TF Recovery time (minutes)", min_value=0, step=1, value=20,
                                key='TF_recovery_time')

    protocol = False
    if st.button("Make Protocol"):
        export_JSON = {
            "Plates": {},
            "Reactions": {},
            "Reaction_volume": {},
            "Deck": {},
            "Parameters": {},
        }
        plate_types = ["source", "Reaction", "TF"]
        for plate_type in plate_types:
            for num in range(st.session_state[f"num_of_{plate_type}_plate"]):
                if st.session_state[f"{plate_type}_plate_{num}_toggle"]:
                    editor_update(
                        editor_key=f"{plate_type}_plate_{num}_long",
                        df=st.session_state[f"{plate_type}_plate_{num}_df"],
                        df_form="long",
                    )
                else:
                    editor_update(
                        editor_key=f"{plate_type}_plate_{num}_wide",
                        df=st.session_state[f"{plate_type}_plate_{num}_df"],
                        df_form="wide",
                    )
                # Json append With Long form
                export_JSON["Plates"][f"{plate_type}_{str(num)}"] = {
                    "name": st.session_state[f"{plate_type}_plate_{num}_name"],
                    "type": plate_type,
                    "labware": st.session_state[f"{plate_type}_plate_{num}_labware"],
                    "data": (
                        plate_transformation(
                            st.session_state[f"{plate_type}_plate_{num}_df"], "long"
                            )
                        .dropna()
                        .reset_index()
                        .set_index("well")).to_dict()["Value"],
                }

        reaction_types = ["PCR", "Assembly"]
        for n, reaction_type in enumerate(reaction_types):
            editor_update(
                editor_key=f"{reaction_type}_plate_{0}_wide",
                df=st.session_state[f"{reaction_type}_plate_{0}_df"],
                df_form="wide",
            )
            export_JSON["Reactions"][f"{reaction_type}_{n}"] = {
                "name": f"{reaction_type}_{n}",
                "type": reaction_type,
                "data": st.session_state[f"{reaction_type}_plate_{0}_df"].astype(str).to_dict()
            }

            # Reaction Volume
            editor_update(
                editor_key=f"{reaction_type}_volume",
                df=eval(f"{reaction_type}_volume"),
                df_form="wide",
            )
            export_JSON["Reaction_volume"][f"{reaction_type}"] = (
                eval(f"{reaction_type}_volume").iloc[0, :].to_dict()
            )

        # Check overlap name in DNA and Reaction
        # DNA overlap from export_JSON["Plates"]
        DNA_names, Reaction_names, TF_names = [], [], []
        for plate_info in export_JSON["Plates"].values():
            if plate_info["type"] == "source":
                DNA_names += list(plate_info["data"].values())
            elif plate_info["type"] == "Reaction":
                Reaction_names += list(plate_info["data"].values())
            elif plate_info["type"] == "TF":
                TF_names += list(plate_info["data"].values())

        # Remove "nan"
        DNA_names = [i for i in DNA_names if i != "nan"]
        Reaction_names = [i for i in Reaction_names if i != "nan"]
        TF_names = [i for i in TF_names if i != "nan"]

        if len(DNA_names) + len(Reaction_names) != len(set(DNA_names + Reaction_names)):
            if len(DNA_names) != len(set(DNA_names)):
                for i in DNA_names:
                    if DNA_names.count(i) > 1:
                        st.error(f"Check DNA names: {i} is Overlaped")
                        
            elif len(Reaction_names) != len(set(Reaction_names)):
                for i in Reaction_names:
                    if Reaction_names.count(i) > 1:
                        st.error(f"Check Reaction names: {i} is Overlaped")
    
            else:
                for i in set(DNA_names) & set(Reaction_names):
                    st.error(f"Check DNA and Reaction names: {i} is Overlaped")
            st.stop()
            
        for i in TF_names:
            if pd.isna(i) or i == "":
                continue
            assert i in Reaction_names or i in DNA_names, f"Check TF names: {i} is not in Reaction Name"
        del TF_names, Reaction_names

        # All DNAs in PCR and Assembly
        PCR_DNA_columns = [
            column_name
            for column_name in st.session_state["PCR_plate_0_df"].columns
            if re.findall(r"DNA(\d+)", column_name)
        ]
        PCR_DNA = st.session_state["PCR_plate_0_df"][PCR_DNA_columns].values
        Assembly_DNA_columns = [
            column_name
            for column_name in st.session_state["Assembly_plate_0_df"].columns
            if re.findall(r"DNA(\d+)", column_name)
        ]
        Assembly_DNA = st.session_state["Assembly_plate_0_df"][Assembly_DNA_columns].values

        # Make list from nested nd.array
        PCR_DNA = [i for j in PCR_DNA for i in j]
        Assembly_DNA = [i for j in Assembly_DNA for i in j]

        for i in list(set(PCR_DNA)):
            if pd.isna(i) or i == "None" or i == "":
                continue
            assert i in DNA_names, f"Check PCR table: {i} is not in DNA plate"
        for i in list(set(Assembly_DNA)):
            if pd.isna(i) or i == "None" or i == "":
                continue
            if i in st.session_state["PCR_plate_0_df"]["Name"].dropna().unique():
                pass
            else:
                assert i in DNA_names, f"Check Assembly table: {i} is not in DNA plate"
        del PCR_DNA_columns, Assembly_DNA_columns, PCR_DNA, Assembly_DNA, DNA_names

        # All Enzymes in PCR and Assembly
        PCR_enzyme_columns = [
            column_name
            for column_name in st.session_state["PCR_plate_0_df"].columns
            if re.findall(r"Enzyme(\d+)", column_name)
        ]
        Assembly_enzyme_columns = [
            column_name
            for column_name in st.session_state["Assembly_plate_0_df"].columns
            if re.findall(r"Enzyme(\d+)", column_name)
        ]
        PCR_enzyme = st.session_state["PCR_plate_0_df"][PCR_enzyme_columns].dropna().values
        Assembly_enzyme = (
            st.session_state["Assembly_plate_0_df"][Assembly_enzyme_columns].dropna().values
        )

        PCR_enzyme = [i for j in PCR_enzyme for i in j]
        Assembly_enzyme = [i for j in Assembly_enzyme for i in j]

        # Protocol 별 검사 (Tip 수, Plate 수, Volume 최대 등)
        export_JSON["Parameters"] = {
            "protocol": "OT-2 cloning",
            "Plate_type": plate_types,
            "Reaction_type": reaction_types,
            "Messenger": st.session_state['messenger'],
            "Stop_between_reactions": st.session_state['stop_between_reactions'],
            "PCR_extension_time": st.session_state['PCR_extension_time'],
            "TF_recovery_time": st.session_state['TF_recovery_time'],
            "number_of_tips": check_tips()
            }

        export_JSON["Deck"] = {
            "Enzyme_position": enzyme_position(
                enzyme_list=list(set(PCR_enzyme) | set(Assembly_enzyme))
            ),
            "Deck_position": deck_position(export_JSON["Plates"]),
        }

        protocol = True
        with open("data/ot2/protocols/cloning/ot2_cloning.py", "r") as f:
            protocol = f.read()
            protocol = protocol.replace("{{EXPORT_JSON}}", str(export_JSON))
            protocol = protocol.replace("{{PRESENT_TIME}}", datetime.now().strftime('%Y-%m-%d'))
            protocol = protocol.replace("nan", '')
            protocol = protocol.replace("None", '')
            new_protocol = protocol.replace("null", '')
            new_protocol=protocol

    st.download_button(
        label="Download Protocol",
        data=new_protocol if protocol else "",
        file_name=f"{datetime.now().strftime('%Y%m%d')[2:]}-ot2_cloning.py",
        disabled=not bool(protocol),
    )
    if export_JSON:
        with st.expander("Report", expanded=True):
            with st.container(height=450):
                st.json(export_JSON)

if __name__ == "__main__":
    st.set_page_config(page_title="Cloning", layout="wide")
    main()
