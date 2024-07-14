import streamlit as st
from streamlit import session_state as state
import pandas as pd
import json
from pathlib import Path

# def
def main():
    def check_project():
        return list(Path("data/project").glob('*.json'))

    def load_project(value: Path):
        with open(value, 'r') as f:
            project_json = json.load(f)
        
        # Duplicate check
        if str(value.stem) in state.loaded_project['Project'].values:
            return
        
        state.loaded_project = \
            pd.concat([state.loaded_project, pd.DataFrame({'Project':[str(value.stem)], 'Task': ', '.join(project_json['Build'].keys()) + ', Transformation'})], axis=0)

    def add_workflow(value):
        tmp = [i.startswith(value) for i in state.new_workflow]
        if not len(tmp):
            num = 1
        else:
            num = len(tmp) + 1
        
        state.new_workflow.append(f'{value}_{num}')

    def reset_workflow():
        state.new_workflow = state.select_workflow

    def state_initiation(key, value):
        if key not in state:
            state[key] = value

    def empty_plate_df(plate_type='96well'):
        if plate_type == '96well':
            rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
            columns = [str(i) for i in range(1, 13)]     
            
            long_columns = ['Value']
            long_rows = [i+j for j in columns for i in rows]
            
            plate_df = pd.DataFrame(index=long_rows, columns=long_columns)
            plate_df.index.name = 'well'
            return plate_df

    def toggle_change(key):
        state_toggle = state[f'{key}_toggle']
        state_edit = state[f'{key}_edit_plate']
        
        if state_toggle:
            state[f'{key}_plate'] = plate_transformation(state_edit, 'wide')
        else:
            state[f'{key}_plate'] = plate_transformation(state_edit, 'long')

    def plate_table(plate_type, use_name=True, loaded_table=False):
        with st.expander(f"{plate_type}", expanded=True):
            st.number_input(f"Number of {plate_type} plate",
                            min_value=1, step=1, max_value=3,
                            key=f"{plate_type}_num")
            
            plates = st.tabs([f"{plate_type}_Plate_{i+1}" for i in range(state[f"{plate_type}_num"])])
            
            for n in range(len(plates)):
                if n > 0:
                    state_initiation(f'{plate_type}_{n+1}_plate', empty_plate_df())
                elif type(loaded_table) != bool:
                    state_initiation(f'{plate_type}_{n+1}_plate', loaded_table)
                else:
                    state_initiation(f'{plate_type}_{n+1}_plate', empty_plate_df())
                state_initiation(f'{plate_type}_{n+1}_edit_plate', None)
                state_initiation(f'{plate_type}_{n+1}_toggle', False)
                
                with plates[n]:
                    if use_name:
                        st.selectbox("Name",
                                     options=[f"{plate_type}_plate_{n+1}", "Eco Part"],
                                     key=f"{plate_type}_plate_{n+1}_name")
                    if st.toggle("Wide form", value=False, key=f'{plate_type}_{n+1}_toggle',
                                on_change=toggle_change,
                                kwargs={'key':f'{plate_type}_{n+1}'}):
                        state[f'{plate_type}_{n+1}_edit_plate'] = st.data_editor(state[f'{plate_type}_{n+1}_plate'], key=f"{plate_type}_plate_{n+1}_editor",
                                                                                num_rows='fixed')
                    else:
                        state[f'{plate_type}_{n+1}_edit_plate'] = st.data_editor(state[f'{plate_type}_{n+1}_plate'], key=f"{plate_type}_plate_{n+1}_editor",
                                                                                num_rows='fixed')

    def workflow_table(workflow, loaded_table=False):
        with st.expander(f'{workflow}', expanded=True):        
            if type(loaded_table) != bool:
                state_initiation(key=f'{workflow}_table', value=loaded_table)
            else:
                state_initiation(key=f'{workflow}_table', value=initial_tables[f'{workflow.split("_")[0]}'])
            
            state[f'{workflow}_edit_table'] = st.data_editor(state[f'{workflow}_table'], key=f'{workflow}_editor',
                        use_container_width=True,
                        hide_index=True,
                        num_rows='dynamic')
            
            if st.button("Add Column", key=f"{workflow}_plate_addcolumn"):
                # state[f'{workflow}_table'] = state[f'{workflow}_edit_table']
                # dataframe add column
                df = state[f'{workflow}_edit_table']
                assert len(df.columns) < 10, "ERROR1: Not allowed more than 10 columns"
                
                df.loc[:, str(len(df.columns)-3)] = [None for _ in range(len(df.index))]
                df.set_index('Name', inplace=True)
                df = df[sorted(df.columns)]
                df.reset_index(inplace=True)
                state[f'{workflow}_table'] = df                
                st.rerun()

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

    def enzyme_position(enzyme_list):
        # 24-well position
        enzyme_list = enzyme_list
        well_position = ["A1", "A2", "A3", "A4", "A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6", "C1", "C2", "C3", "C4", "C5", "C6", "D1", "D2", "D3", "D4", "D5", "D6"]
        return_dict = {}
        for enzyme, well in zip(enzyme_list, well_position):
            return_dict[enzyme] = well
        return return_dict

    def deck_position(plates, additional_plate: list):
        position = [1,2,3,4,5,6,9]
        
        deck_dict = {}
        deck_dict["Enzyme_tube"] = position.pop(0)
        deck_dict["p20_tip"] = position.pop(0)
        deck_dict["p300_tip"] = position.pop(0)
        
        for key in plates.keys():
            # OT-2
            if plates[key]["type"] == "Destination":
                deck_dict[key] = 7
                continue
            try:
                deck_dict[key] = position.pop(0)
            except:
                raise "Deck is already Full. Reduce Plates"
        
        if additional_plate:
            for key in additional_plate:
                try:
                    deck_dict[key] = position.pop(0)
                except:
                    raise "Deck is already Full. Reduce Plates"
        
        return deck_dict



    # Statics
    if 'project' not in state:
        state.project = check_project()
    if 'new_workflow' not in state:
        state.new_workflow = []
    if 'loaded_project' not in state:
        state.loaded_project = pd.DataFrame({'Project':[], 'Task':[]}, index=None)
    if 'make_workflows' not in state:
        state.make_workflows = False
    if 'export_JSON' not in state:
        state.export_JSON = False
    if 'make_json' not in state:
        state.make_json = False

    workflows = ['PCR', 'Gibson', 'GGA', 'TF']
    initial_tables = {
        "empty": empty_plate_df(),
        "PCR": pd.DataFrame({
                    "Name": [None],
                    "0": [None],
                    "1": [None],
                    "2": [None],
                    "A_enzyme": ["[E]PCRmix"],
                    "DW": ["[E]DW"],
                }),
        "Gibson": pd.DataFrame({
                    "Name": [None],
                    "0": [None],
                    "1": [None],
                    "A_enzyme": ["[E]Gibsonmix"],
                    "DW": ["[E]DW"],
                }),
        "GGA": pd.DataFrame({
                    "Name": [None],
                    "0": [None],
                    "1": [None],
                    "2": [None],
                    "3": ["[E]BsaI"],
                    "4": ["[E]T4_ligase"],
                    "A_enzyme": ["[E]Buffer"],
                    "DW": ["[E]DW"],
                })
        }

    # Main
    col1 = st.columns([1, 2])
    with col1[0]:
        st.markdown("## Load Project")
        st.selectbox("Saved Project", state.project,
                    format_func=lambda x: x.stem,
                    key='select_project')

        st.button("Load", key='load_project',
                type='secondary',
                on_click=load_project,
                kwargs={'value': state.select_project})

    no_use_new_task = bool(len(state.loaded_project))
    with col1[1]:
        st.markdown("## New Task")
        st.multiselect("New workflows", state['new_workflow'],
                    key='select_workflow',
                    default=state.new_workflow,
                    on_change=reset_workflow,
                    help='New Task is dominant than load project')
        workflow_col = st.columns([1,1,1,1,1])
        workflow_col[0].button("PCR", on_click=add_workflow,
                            kwargs={'value': 'PCR'},
                            use_container_width=True, disabled=no_use_new_task)
        workflow_col[1].button("GGA", use_container_width=True,
                            on_click=add_workflow,
                            kwargs={'value': 'GGA'}, disabled=no_use_new_task)
        workflow_col[2].button("Gibson", use_container_width=True,
                            on_click=add_workflow,
                            kwargs={'value': 'Gibson'}, disabled=no_use_new_task)
        workflow_col[3].button("TF", use_container_width=True,
                            on_click=add_workflow,
                            kwargs={'value': 'Transformation'}, disabled=no_use_new_task)
        
        if state.select_workflow or len(state.loaded_project.index):
            use_make = True
        else:
            use_make = False
        
        if workflow_col[4].button("Make", use_container_width=True, 
                                type="primary",
                                disabled=not use_make,
                                help="When Project loaded or Make new workflows open."):
            if state.make_workflows:
                st.error("Workflows are already loaded please Restart App.")
            state.make_workflows = True

    st.file_uploader("## Load previous result", type=['py'],key="load", disabled=True, help="Load previous result to modify")
    st.dataframe(state.loaded_project, use_container_width=True, hide_index=True)
    st.markdown('---')

    mid_col = st.columns([1, 1])
    with mid_col[0]:
        st.markdown('## Mapping')
        
        # 호출 시 Project 데이터를 받아오기 위해 필요함
        if state.make_workflows:
            source, dest = [], []
            for project in state.loaded_project['Project']:
                with open(f'data/project/{project}.json', 'r') as f:
                    js = json.load(f)
                    for key in js['Build'].keys():
                        source += list(js['Build'][key].values())
                        dest += list(js['Build'][key].keys())
            try:
                source = sum(source, [])
                dest = sum(dest, [])
            except:
                pass
            source = list(dict.fromkeys(source))
            dest = list(dict.fromkeys(dest))
            
            source = [i for i in source if not i.startswith('[E]')]
            source = [i for i in source if i not in dest]
            
            source_df = empty_plate_df()
            source_df.iloc[:len(source), 0] = source
            dest_df = empty_plate_df()
            # 짝수 줄 건너뛰게 해야 함
            dest_df.iloc[:len(dest), 0] = dest
            
            # Source Plate Module
            plate_table("Source", use_name=True, loaded_table=source_df)
            plate_table("Destination", use_name=False, loaded_table=dest_df)

    with mid_col[1]:
        st.markdown('## Workflow')
        if state.make_workflows:
            # no loaded project
            if not no_use_new_task:
                state.workflow = state.select_workflow
                for workflow in state.workflow:
                    if workflow.split('_')[0] == 'Transformation':
                        plate_table(workflow, use_name=True)
                    else:
                        workflow_table(workflow) # 여기서는 정상
            # Use loaded Project
            else:
                # 로드된 것 중 가장 긴 프로젝트 선택
                tmp = state.loaded_project['Task'].apply(lambda x:x.split(', '))
                cnt = 0
                for i in tmp:
                    if len(i) > cnt:
                        cnt = len(i)
                        state.workflow = i
                
                # Project에 처음부터 Workflow_1 으로 기록해야 할 듯..
                state.workflow=[i+f'_{str(n)}' for n, i in enumerate(state.workflow, 1)]
                
                for workflow in state.workflow:
                    if workflow.split('_')[0] == 'Transformation':
                        plate_table(workflow, use_name=True)
                    else:
                        # Load table by workflow
                        df = pd.DataFrame()
                        for project in state.loaded_project['Project']:
                            with open(f'data/project/{project}.json', 'r') as f:
                                js = json.load(f)
                                if not workflow.split('_')[0] in js['Build'].keys():
                                    continue
                                new_df = pd.DataFrame(js['Build'][workflow.split('_')[0]]).T.reset_index(names=["Name"])
                                df = pd.concat([df, new_df], axis=0, ignore_index=True)
                        workflow_table(workflow, loaded_table=df)

    # Parameters
    st.markdown('---')
    st.markdown('## Parameters')
    with st.expander("Parameters", expanded=True):
        if state.make_workflows:
            st.success("Please adjust here as last step")
            # Reaction-PCR: DNA, Enzyme, DW(up to)
            # Reaction-Assembly: DNA, Enzyme, DW(up to)

            for workflow in state.workflow:
                if workflow.startswith('Transformation'):
                    continue
                st.markdown(f'### {workflow} volume')
                df = (state[f"{workflow}_edit_table"]
                    .drop(["Name"], axis=1, inplace=False)
                    .dropna(axis=1, how='all')
                    .reset_index(drop=True, inplace=False))
                
                # remove all values in df
                df = df.iloc[[0], :]
                df.loc[:, :] = None
                
                if workflow.startswith('PCR'):
                    df.loc[0, ['0', '1', '2', 'A_enzyme', 'DW']] = ['1', '0.5', '0.5', '12.5', '10.5']
                elif workflow.startswith('Gibson'):
                    df.loc[0, ['0', '1', 'A_enzyme', 'DW']] = ['2', '2', '5', '1']
                elif workflow.startswith('GGA'):
                    df.loc[0, ['0', '1', '2', '3', '4', 'A_enzyme', 'DW']] = ['1', '1', '1', '1', '0.5', '2.5', '3']
                
                df.index.name = 'Index'
                df = df[sorted(df.columns)]
                state[f"{workflow}_edit_volume"] = st.data_editor(df, use_container_width=True, key=f"{workflow}_volume",
                                                                  hide_index=True,
                                                                  num_rows='fixed')
            
            
            st.checkbox('Stop between Reactions', value=True,
                        key='stop_reaction',
                        help='If you want to stop between reactions for take enzymes, check this box')

            advanced_column = st.columns([1,1])
            with advanced_column[0]:
                with st.container(border=True):
                    st.number_input("Annealing temperature", min_value=45, step=1, value=57,
                                    key='annealing')
                    st.number_input("PCR extension time (seconds)", min_value=1, step=1, value=25,
                                    key='pcr_extension')
                    st.number_input("TF Recovery time (minutes)", min_value=0, step=1, value=40,
                                    key='tf_recovery')

    end_col = st.columns([1,1])
    with end_col[0]:
        if st.button("Make Protocol", type="primary",
                    disabled = not state.make_workflows):
            state.make_json = True
        
        if state.make_json:
            state.export_JSON = {
                "Meta": {},
                "Plate": {},
                "Workflow": {},
                "Workflow_volume": {},
                "Deck": {},
                "Parameter": {},
            }
            
            # Meta
            state.export_JSON["Meta"] = {
                "Task": "OT-2 cloning",
                "version": "2.1",
                "workflow": state.workflow,
                "Messenger": "kun"
            }
            plate_types = ["Source", "Destination"]
            for plate_type in plate_types:
                for n in range(state[f"{plate_type}_num"]):
                    if state[f'{plate_type}_{n+1}_toggle']:
                        # Wide to Long
                        value = plate_transformation(state[f'{plate_type}_{n+1}_edit_plate'], 'long').dropna()["Value"].to_dict()
                    else:
                        value = state[f'{plate_type}_{n+1}_edit_plate'].dropna()["Value"].to_dict()
                    if f"{plate_type}_plate_{n+1}_name" in state:
                        name = state[f"{plate_type}_plate_{n+1}_name"]
                    else:
                        name = f"{plate_type}_{n+1}"
            
                    state.export_JSON["Plate"][f"{plate_type}_{n+1}"] = {
                        "name": name,
                        "type": plate_type,
                        "data": value,
                    }
            
            # Workflow
            for workflow in state.workflow:
                if workflow.startswith("Transformation"):
                    for n in range(state[f"{workflow}_num"]):
                        if f"{workflow}_plate_{n+1}_name" in state:
                            name = state[f"{workflow}_plate_{n+1}_name"]
                        else:
                            name = f"{workflow}_{n+1}"
                        if state[f'{workflow}_{n+1}_toggle']:
                            value = plate_transformation(state[f'{workflow}_{n+1}_edit_plate'], 'long').dropna()["Value"].to_dict()
                        else:
                            value = state[f'{workflow}_{n+1}_edit_plate'].dropna()["Value"].to_dict()
                        # workflow가 여러개가 들어가는 형태로 되었음.. data가 여러개가 들어가야 할 것 같은뎅
                        state.export_JSON["Workflow"][f"{workflow}_{n+1}"] = {
                            "name": name,
                            "type": workflow.split('_')[0],
                            "data": value,
                        }
                else:
                    value = state[f'{workflow}_edit_table'].astype(str).to_dict()

                    # Streamlit 자체 이슈로 변환 과정 중 sort가 걸림.
                    state.export_JSON['Workflow'][workflow] = {
                        "type": workflow.split('_')[0],
                        "data": value
                    }
            
            # Workflow volume
            for workflow in state.workflow:
                if workflow.startswith("Transformation"):
                    continue
                
                name = workflow
                value = state[f"{workflow}_edit_volume"].astype(str).to_dict()
                # None 이 있으면 Error 발생
                assert None not in state[f"{workflow}_edit_volume"].values, "ERROR2: Fill, all of Volume tables!"
                assert "" not in state[f"{workflow}_edit_volume"].values, "ERROR2: Fill, all of Volume tables!"
            
                state.export_JSON['Workflow_volume'][name] = value
            
            # Parameter
            state.export_JSON["Parameter"] = {
                "stop_reaction": state.stop_reaction,
                "annealing": state.annealing,
                "pcr_extension": state.pcr_extension,
                "tf_recovery": state.tf_recovery,
                "num_of_tips": "NULL"
            }
            
            # Deck
            ## Materials
            use_tf = False
            materials, enzymes, dnas, products = [], [], [], []
            for workflow in state.workflow:
                if workflow.startswith("Transformation"):
                    use_tf = True
                    continue
                tmp = state[f'{workflow}_edit_table']
                products += tmp["Name"].tolist()
                materials += tmp.drop(["Name"], axis=1).values.tolist()
            try:
                materials = sum(materials, [])
                products = sum(products, [])
            except:
                pass
            
            for i in list(dict.fromkeys(materials)):
                if type(i) != str:
                    continue
                if i.startswith('[E]'):
                    enzymes.append(i)
                else:
                    dnas.append(i)
            
            if use_tf:
                enzymes += ["CPcell", "SOC"]
            
            ## Deck position
            tf_plate = []
            for key in state.export_JSON["Workflow"].keys():
                if state.export_JSON["Workflow"][key]["type"] == "Transformation":
                    tf_plate.append(key)
            
            state.export_JSON["Deck"] = {
                "Enzyme_position": enzyme_position(enzymes),
                "Deck_position": deck_position(state.export_JSON["Plate"], tf_plate)
            }
            
            # Check error
            source_materials = []
            dest_materials = []
            
            for key in state.export_JSON["Plate"].keys():
                tmp = state.export_JSON["Plate"][key]
                
                if tmp["type"] == "Source":
                    source_materials += list(tmp["data"].values())
                elif tmp["type"] == "Destination":
                    dest_materials += list(tmp["data"].values())

            try:
                source_materials.remove("")
                dest_materials.remove("")
            except:
                pass
            
            ## DEST에 Product 검사 
            # Product 중복 검사
            assert len(set(products)) == len(products), "ERROR3: Duplicated Product Exists! or Empty"
            # DEST table에 REACTION - NAME들 있는지.
            for element in products:
                if element == None:
                    continue
                assert element in dest_materials, f"ERROR4: Product `{element}` not in Destination Plate!"
            
            ## Source 검사
            # Source or dest plate에 DNA들 있는지 확인
            for element in dnas:
                if element == "":
                    continue
                assert (element in source_materials) or (element in dest_materials), f"ERROR5: source `{element}` not in Plate"
            # Source DNA가 Product와 이름이 겹칠 때
            
            ## TF 검사
            # TF 대상이 Product or Source에 없을 때
            tf_products = []
            for key in tf_plate:
                tf_products += list(state.export_JSON["Workflow"][key]['data'].values())
            
            for element in tf_products:
                assert element in dest_materials, f"ERROR6: TF_product `{element}` not in Destination Plate!"
            
            ## Overlap
            ## 같은 이름이 존재할 때
            ## TF 대상 Product에 없을 때
            state.make_json = False

    if state.export_JSON:
        with st.expander("Converted JSON", expanded=True):
            with st.container(height=450):
                st.json(state.export_JSON)

    with end_col[1]:
        st.download_button(
            label = "Download Protocol",
            data = "",
            disabled=True,
            help="Blocked by DEMO version"
        )