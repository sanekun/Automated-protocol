from opentrons import protocol_api
from opentrons import simulate

protocol = simulate.get_protocol_api("2.13")
import pandas as pd
import re
import time
import json
import requests
import logging

metadata = {
    "protocolName": "{{PRESENT_TIME}} Cloning (PCR, Assembly, Transformation)",
    "author": "Seong-Kun Bak <sanekun@kaist.ac.kr>",
    "apiLevel": "2.13",
    "robotType": "OT-2",
    "description": "Cloning Protocol in SBL",
}

PARAMETERS = {
    "Meta": {
        "Task": "OT-2 cloning",
        "version": "2.1",
        "workflow": ["PCR_1", "GGA_2", "Gibson_3", "Transformation_4"],
        "Messenger": "kun",
    },
    "Plate": {
        "Source_1": {
            "name": "Source_plate_1",
            "type": "Source",
            "data": {"A1": "a", "B1": "s", "C1": "d"},
        },
        "Destination_1": {
            "name": "Destination_1",
            "type": "Destination",
            "data": {"A1": "asdf", "B1": "bde", "C1": "final"},
        },
        "Transformation_4_1": {
            "name": "240715_Transformation_4_1",
            "type": "Transformation",
            "data": {"A1": "final"},
        },
        "Transformation_4_2": {
            "name": "Cm",
            "type": "Transformation",
            "data": {"A1": "bde"},
        },
    },
    "Workflow": {
        "PCR_1": {
            "type": "PCR",
            "data": {
                "0": {"0": "a"},
                "1": {"0": "s"},
                "2": {"0": "d"},
                "Name": {"0": "asdf"},
                "A_enzyme": {"0": "[E]PCRmix"},
                "DW": {"0": "[E]DW"},
            },
        },
        "GGA_2": {
            "type": "GGA",
            "data": {
                "0": {"0": "asdf"},
                "1": {"0": "a"},
                "2": {"0": "d"},
                "3": {"0": "[E]BsaI"},
                "4": {"0": "[E]T4_ligase"},
                "Name": {"0": "bde"},
                "A_enzyme": {"0": "[E]Buffer"},
                "DW": {"0": "[E]DW"},
            },
        },
        "Gibson_3": {
            "type": "Gibson",
            "data": {
                "0": {"0": "bde"},
                "1": {"0": "asdf"},
                "Name": {"0": "final"},
                "A_enzyme": {"0": "[E]Gibsonmix"},
                "DW": {"0": "[E]DW"},
            },
        },
    },
    "Workflow_volume": {
        "PCR_1": {
            "0": {"0": "1"},
            "1": {"0": "0.5"},
            "2": {"0": "0.5"},
            "A_enzyme": {"0": "12.5"},
            "DW": {"0": "10.5"},
        },
        "GGA_2": {
            "0": {"0": "1"},
            "1": {"0": "1"},
            "2": {"0": "1"},
            "3": {"0": "1"},
            "4": {"0": "0.5"},
            "A_enzyme": {"0": "2.5"},
            "DW": {"0": "3"},
        },
        "Gibson_3": {
            "0": {"0": "2"},
            "1": {"0": "2"},
            "A_enzyme": {"0": "5"},
            "DW": {"0": "1"},
        },
    },
    "Deck": {
        "Enzyme_position": {
            "[E]PCRmix": "A1",
            "[E]DW": "A2",
            "[E]BsaI": "A3",
            "[E]T4_ligase": "A4",
            "[E]Buffer": "A5",
            "[E]Gibsonmix": "A6",
            "[E]CPcell": "B1",
            "[E]SOC": "B2",
        },
        "Deck_position": {
            "Enzyme_tube": 1,
            "p20_tip": 2,
            "p300_tip": 3,
            "Source_1": 4,
            "Destination_1": 7,
            "Transformation_4_1": 5,
            "Transformation_4_2": 6,
        },
    },
    "Parameter": {
        "stop_reaction": True,
        "annealing": 57,
        "pcr_extension": 25,
        "tf_recovery": 40,
        "num_of_tips": "NULL",
    },
}

default_labware = "biorad_96_wellplate_200ul_pcr"
logging.basicConfig(filename=f"{metadata['protocolName']}.log", level=logging.INFO)
logging.info(f"Protocol Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
logging.info(f"user: {PARAMETERS['Meta']['Messenger']}")

def run(protocol: protocol_api.ProtocolContext):

    # [Functions]
    def discord_message(message, user=PARAMETERS["Meta"]["Messenger"]):
        if user == "None":
            return None

        # Send message to discord chaneel
        url = "https://discord.com/api/webhooks/1169608986891390996/wLti5ulSOXBtGelwiYO4SJi1jQSJ9tEQ8uC8sVXVlWXih3XWwkoLDr6cJBEm2iPY9b0t"

        headers = {"Content-Type": "application/json"}
        data = {"content": message}
        response = requests.post(url, headers=headers, data=json.dumps(data), verify=False)


    def flow_rate(pipette, **kwargs):
        # Change flow rate of pipette

        assert (
            item in ["aspirate", "dispense", "blow_out"] for item in kwargs.keys()
        ), "Error Keywords in Flow Rate."
        for i in kwargs.keys():
            setattr(pipette.flow_rate, i, kwargs[i])


    def find_materials_well(material, right_well=False):
        # Convert material name to well
        # right_well means well of right side of plate (for abstraction)
        if material.startswith("[E]"):
            return PARAMETERS["Deck"]["Enzyme_position"][material]

        for plate in PARAMETERS["Plate"].values():
            for well in plate["data"].keys():
                if plate["data"][well] == material:
                    wells = plate["Deck"].wells_by_name()
                    if right_well:
                        pos = list(wells.keys()).index(well)
                        return wells[list(wells)[pos + 8]]
                    else:
                        return wells[well]


    def transfer_materials(workflow_df, volume_dict, mix_last=(0, 0)):
        # transform dict to dataframe (row iterable)
        df = pd.DataFrame(workflow_df["data"])
        # Only First Material transfer to all wells at once (Enzyme1 or DW)
        for dw in df["DW"].unique():
            tmp = df[df["DW"] == dw]
            # If Data doens't exist, this step will be skipped
            if len(tmp):
                if pd.isna(dw) or dw == "":
                    continue
                src = find_materials_well(dw)
                vol = float(volume_dict["DW"])
                dest = [find_materials_well(name) for name in tmp["Name"].values]

                # volume에 따라 tip을 달리 사용하도록 하기.
                flow_rate(p300, aspirate=50, dispense=50, blow_out=20)
                p300.distribute(
                    vol,
                    src,
                    dest,
                    new_tip="once",
                    touch_tip=False,
                    disposal_volume=5,
                    blow_out=False,
                )

        for enzyme_name in df["A_enzyme"].unique():
            tmp = df[df["A_enzyme"] == enzyme_name]
            # If Data doens't exist, this step will be skipped
            if len(tmp):
                if pd.isna(enzyme_name) or enzyme_name == "":
                    continue
                src = find_materials_well(enzyme_name).bottom(z=3)
                vol = float(volume_dict["A_enzyme"])
                dest = [find_materials_well(name) for name in tmp["Name"].values]

                flow_rate(p300, aspirate=20, dispense=20, blow_out=20)
                p300.distribute(
                    vol,
                    src,
                    dest,
                    new_tip="once",
                    touch_tip=False,
                    mix_before=(2, 50),
                    disposal_volume=5,
                    blow_out=True,
                    blowout_location="source well",
                )
        # Other Materials
        df.drop(columns=["A_enzyme", "DW"], inplace=True)
        columns = df.columns
        for value in df.values:
            for sample_type, sample_name in zip(columns, value):
                # Empty well will be skipped
                if sample_type == "Name":
                    dest = find_materials_well(sample_name)
                    continue
                if pd.isna(sample_name) or sample_name == "":
                    continue

                src = find_materials_well(sample_name)
                vol = float(volume_dict[sample_type])

                # Check DNA or Enzyme
                if sample_name.startswith("[E]"):
                    flow_rate(p20, aspirate=1, dispense=1, blow_out=1)
                    touch_tip = False
                else:  # DNA & DW
                    flow_rate(p20, aspirate=1, dispense=1, blow_out=1)
                    touch_tip = False

                p20.transfer(
                    vol,
                    src,
                    dest,
                    new_tip="always",
                    touch_tip=touch_tip,
                    blow_out=False,
                    blowout_location="destination well",
                )

            # Mix Product
            if sum(mix_last):
                flow_rate(p20, aspirate=10, dispense=10, blow_out=10)
                p20.pick_up_tip()
                for _ in range(mix_last[0]):
                    p20.aspirate(mix_last[1], dest.bottom())
                    p20.dispense(mix_last[1], dest.bottom(z=3))
                p20.drop_tip()


    def spotting_dispense(pipette, src, dest: list, spotting_volume=4):
        # Dispense liquid at the top of well
        # and then, move down pipette to specific height.

        if not pipette.has_tip:
            pipette.pick_up_tip()

        disposal_vol = 1
        whole_vol = spotting_volume * len(dest) + disposal_vol
        cnt = 0
        while whole_vol > spotting_volume:
            if pipette.max_volume < whole_vol:
                pipette.aspirate(pipette.max_volume, src)
            else:
                pipette.aspirate(whole_vol, src)

            while pipette.current_volume > spotting_volume:
                pipette.dispense(spotting_volume, dest[cnt].bottom(z=4.4))
                pipette.move_to(dest[cnt].bottom(z=3))
                cnt += 1
                whole_vol -= spotting_volume

            pipette.blow_out(pipette.trash_container.wells()[0])


    def run_PCR(workflow_df, volume_dict):
        final_volume = sum(map(float, volume_dict.values()))
        transfer_materials(workflow_df=workflow_df, volume_dict=volume_dict, mix_last=(2, 15))
        ## Thermocycler
        discord_message(f"Thermocycler in {workflow} start RUN take off Enzyme")
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(95)
        tc_mod.set_block_temperature(
            temperature=94, hold_time_seconds=30, block_max_volume=final_volume
        )
        profile = [
            {"temperature": 94, "hold_time_seconds": 20},
            {
                "temperature": int(PARAMETERS["Parameter"]["annealing"]),
                "hold_time_seconds": 20,
            },
            {
                "temperature": 68,
                "hold_time_seconds": int(PARAMETERS["Parameter"]["pcr_extension"]),
            },
        ]
        tc_mod.execute_profile(steps=profile, repetitions=30, block_max_volume=final_volume)
        tc_mod.set_block_temperature(
            temperature=68, hold_time_seconds=60, block_max_volume=final_volume
        )
        discord_message(f"{workflow} will be end 5 minutes later, Take in Enzyme for next step")
        tc_mod.set_block_temperature(
            temperature=12, hold_time_minutes=5, block_max_volume=final_volume
        )
        tc_mod.deactivate_lid()
        tc_mod.set_block_temperature(12)
        # End with closed lid

    def run_Gibson(workflow_df, volume_dict):
        final_volume = sum(map(float, volume_dict.values()))

        transfer_materials(workflow_df=workflow_df, volume_dict=volume_dict, mix_last=(2, 15))
        discord_message(f"{key}: Thermocycler is running remove Enzyme")

        ## Thermocycler
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(80)

        ### DpnI
        tc_mod.set_block_temperature(
            temperature=37, hold_time_minutes=5, block_max_volume=final_volume
        )
        ### denaturation
        tc_mod.set_block_temperature(
            temperature=65, hold_time_seconds=20, block_max_volume=final_volume
        )
        tc_mod.set_block_temperature(
            temperature=50,
            hold_time_minutes=40,
            block_max_volume=final_volume,
        )
        discord_message(
            f"{workflow} will be end 10 minutes later, Take in CP cell for next step"
        )
        # Ramp rate is 0.1 degree per second
        current_tmp = 50
        while True:
            current_tmp -= 5

            if current_tmp <= 12:
                tc_mod.deactivate_lid()
                tc_mod.set_block_temperature(12)
                break

            tc_mod.set_block_temperature(
                temperature=current_tmp,
                hold_time_seconds=45,  # ramp_rate is almost 0.1 degree per second
                block_max_volume=final_volume,
            )


    def run_GGA(workflow_df, volume_dict):
        final_volume = sum(map(float, volume_dict.values()))
        transfer_materials(workflow_df=workflow_df, volume_dict=volume_dict, mix_last=(2, 15))
        ## Thermocycler
        discord_message(f"{key}: Thermocycler in PCR start RUN take off Enzyme")

        tc_mod.close_lid()
        tc_mod.set_lid_temperature(90)

        profile = [
            {"temperature": 37, "hold_time_seconds": 20},
            {"temperature": 16, "hold_time_seconds": 20},
        ]
        tc_mod.execute_profile(steps=profile, repetitions=30, block_max_volume=final_volume)

        discord_message(f"{workflow} will be end 5 minutes later, Take in Enzyme for next step")
        tc_mod.set_block_temperature(
            temperature=12, hold_time_minutes=5, block_max_volume=final_volume
        )
        tc_mod.deactivate_lid()
        tc_mod.set_block_temperature(12)


    def run_Transformation():
        flow_rate(p300, aspirate=20, dispense=20, blow_out=100)
        src = find_materials_well("[E]CPcell")

        unique_sample = []
        for key in PARAMETERS["Plate"].keys():
            if not key.startswith("Transformation"):
                continue

            plate = PARAMETERS["Plate"][key]
            if plate["type"] == "Transformation":
                unique_sample += list(plate["data"].values())
        unique_sample = list(set(unique_sample))
        if "" in unique_sample:
            unique_sample.remove("")

        dest = [
            # use right well of destination plate for CP cell
            find_materials_well(name, right_well=True)
            for name in unique_sample
        ]

        CP_cell_volume = 45
        ## Mix CP cell
        p300.pick_up_tip()
        for _ in range(2):
            p300.aspirate(25, src)
            p300.dispense(25, src.bottom(z=10))

        # Transfer CP cell to right well of assembly product
        p300.distribute(
            CP_cell_volume,
            src.bottom(z=3),
            dest,
            new_tip="never",
            touch_tip=False,
            disposal_volume=10,
            blow_out=True,
            blowout_location="source well",
        )
        p300.drop_tip()

        # Transfer Assembly Mix to distributed CP cell
        src = [find_materials_well(name, "DNA", right_well=False) for name in unique_sample]
        reaction_mix_vol = 5
        p20.transfer(reaction_mix_vol, src, dest, new_tip="always", blow_out=False)

        tc_mod.close_lid()
        protocol.delay(minutes=10)

        ## Thermocycler
        tc_mod.set_block_temperature(
            temperature=42,
            hold_time_seconds=90,
            block_max_volume=reaction_mix_vol + CP_cell_volume,
        )
        tc_mod.set_block_temperature(8)
        tc_mod.open_lid()

        src = find_materials_well("LB", "DW").bottom(z=3)

        # Add media for recovery
        # start_time = time.time()
        p300.transfer(
            100,
            src,
            dest,
            new_tip="always",
            touch_tip=False,
            disposal_volume=5,
            blow_out=False,
        )
        protocol.delay(seconds=30)
        # end_time = time.time()

        # Duration time in 8 degree
        # rest_time = 120 - (end_time - start_time)
        # if rest_time > 0:
        #    if rest_time < 30:
        #        protocol.delay(seconds=30)
        #    else:
        #        protocol.delay(seconds=rest_time)

        # Recovery
        tc_mod.set_block_temperature(
            temperature=37,
            hold_time_minutes=int(int(PARAMETERS["Parameters"]["TF_recovery_time"]) / 2),
        )

        for dest_well in dest:
            p300.pick_up_tip()
            for _ in range(2):
                p300.aspirate(40, dest_well)
                p300.dispense(40, dest_well.bottom(z=5))
            p300.drop_tip()

        tc_mod.set_block_temperature(
            temperature=37,
            hold_time_minutes=int(int(PARAMETERS["Parameters"]["TF_recovery_time"]) / 2),
        )


        # Spotting
        spotting_volume = 4
        flow_rate(p20, aspirate=8, dispense=15, blow_out=15)
        for i in PARAMETERS["Plates"].keys():
            plate = PARAMETERS["Plates"][i]
            if plate["type"] != "TF":
                continue
            else:
                unique_sample = list(set(plate["data"].values()))
                if "" in unique_sample:
                    unique_sample.remove("")
                for sample in unique_sample:
                    src = find_materials_well(sample, "DNA", right_well=True)
                    dest = [
                        plate["Deck"][well]
                        for well, value in plate["data"].items()
                        if value == sample
                    ]
                    p20.pick_up_tip()
                    # Mix Sample
                    for _ in range(3):
                        p20.aspirate(20, src)
                        p20.dispense(20, src.bottom(z=4))
                    # Dispense Sample
                    spotting_dispense(p20, src, dest, spotting_volume)
                    # p20.distribute(spotting_volume, src, dest, new_tip="never", touch_tip=False)
                    p20.drop_tip()
        tc_mod.deactivate()
    
    #------------------------------------------------ Protocol Start
    discord_message(f"Protocol Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    # Deck Setting
    ## Modules
    tc_mod = protocol.load_module(module_name="thermocyclerModuleV1")
    tc_mod.open_lid()
    # enzyme deck is fixed in ot-2
    Enzyme_deck = protocol.load_labware("opentrons_24_tuberack_nest_1.5ml_screwcap", 1)

    ## Pipette
    p20_tip = protocol.load_labware("opentrons_96_tiprack_20ul", 2)
    p300_tip = protocol.load_labware("opentrons_96_tiprack_300ul", 3)

    p20 = protocol.load_instrument("p20_single_gen2", "left", tip_racks=[p20_tip])
    p300 = protocol.load_instrument("p300_single_gen2", "right", tip_racks=[p300_tip])

    ## Enzymes
    for key in PARAMETERS["Deck"]["Enzyme_position"].keys():
        PARAMETERS["Deck"]["Enzyme_position"][key] = Enzyme_deck.wells_by_name()[
            PARAMETERS["Deck"]["Enzyme_position"][key]
        ]

    ## Plates
    for key in PARAMETERS["Plate"].keys():
        location = PARAMETERS["Deck"]["Deck_position"][key]
        if str(location) == "7":
            PARAMETERS["Plate"][key]["Deck"] = tc_mod.load_labware(default_labware)
            continue

        # TF plate 넣어야 함.
        PARAMETERS["Plate"][key]["Deck"] = protocol.load_labware(
            default_labware, location=location
        )

    ## Workflows
    for workflow in PARAMETERS["Meta"]["workflow"]:
        key = workflow.split('_')[0]
        assert key in [
            "PCR",
            "GGA",
            "Gibson",
            "Transformation",
        ], f"{workflow}: Error Workflow"
        
        # Empty workflow를 무시하고 지나갈 수 있어야 함.
        if PARAMETERS["Parameter"]["stop_reaction"]:
            # 첫 번째 workflow 전은 stop하지 않음
            if not workflow == PARAMETERS["Meta"]["workflow"][0]:
                protocol.pause(f"{workflow}: will be start Place down enzyme")
                discord_message(f"{workflow}: Protocol Paused please push start button")        
        
        if key == "Transformation":
            run_Transformation(workflow_df=PARAMETERS["Workflow"][workflow], volume_dict=PARAMETERS["Workflow_volume"][workflow])
            continue

        workflow_df = PARAMETERS["Workflow"][workflow]
        volume_dict = PARAMETERS["Workflow_volume"][workflow]

        # key - value 형식으로 변경
        for i in volume_dict.keys():
            volume_dict[i] = next(volume_dict[i].values().__iter__())
            
        f"run_{key}(workflow_df={workflow_df},volume_dict={volume_dict})"
        # Run workflow functions
        eval(f"run_{key}({workflow_df},{volume_dict})")
        tc_mod.open_lid()

    discord_message(f"Protocol End: {time.strftime('%Y-%m-%d %H:%M:%S')}")