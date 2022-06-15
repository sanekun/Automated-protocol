from opentrons import types, protocol_api
from math import floor

#|# DNA assembly template script. #|#
#|# Don't Run This script directly with OT2! #|#

metadata = {
    'protocolName': 'Golden gate assembly used by SBL (96 well plate based)',
    'author': 'Seong-Kun Bak <tjdrns27@kribb.re.kr>',
    'apiLevel': '2.11',
    'description': 'Simultaneous part assembly with different part combination',
    'date' : '06/15/22'
}

# Load assembled well data
meta_data = [{'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A1', 'vol': 0.6, 'dil': 0, 'plate': 'pro'}, 'part3': {'well': 'A1', 'vol': 0.55, 'dil': 0, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p1_r1_t25', 'name': 'pACBB_DmpR_BBa_I14018_BBa_J61100_L3S1P56', 'DW': 5.41}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A1', 'vol': 0.6, 'dil': 0, 'plate': 'pro'}, 'part3': {'well': 'B6', 'vol': 2.45, 'dil': 5, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p1_r18_t25', 'name': 'pACBB_DmpR_BBa_I14018_BBa_B0072_L3S1P56', 'DW': 3.51}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A1', 'vol': 0.6, 'dil': 0, 'plate': 'pro'}, 'part3': {'well': 'B10', 'vol': 2.1, 'dil': 5, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p1_r22_t25', 'name': 'pACBB_DmpR_BBa_I14018_BBa_B0032_L3S1P56', 'DW': 3.86}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A1', 'vol': 0.6, 'dil': 0, 'plate': 'pro'}, 'part3': {'well': 'B11', 'vol': 0.52, 'dil': 0, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p1_r23_t25', 'name': 'pACBB_DmpR_BBa_I14018_BBa_B0034_L3S1P56', 'DW': 5.44}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A6', 'vol': 0.53, 'dil': 0, 'plate': 'pro'}, 'part3': {'well': 'A1', 'vol': 0.55, 'dil': 0, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p6_r1_t25', 'name': 'pACBB_DmpR_BBa_J23112_BBa_J61100_L3S1P56', 'DW': 5.48}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A6', 'vol': 0.53, 'dil': 0, 'plate': 'pro'}, 'part3': {'well': 'B6', 'vol': 2.45, 'dil': 5, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p6_r18_t25', 'name': 'pACBB_DmpR_BBa_J23112_BBa_B0072_L3S1P56', 'DW': 3.58}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A6', 'vol': 0.53, 'dil': 0, 'plate': 'pro'}, 'part3': {'well': 'B10', 'vol': 2.1, 'dil': 5, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p6_r22_t25', 'name': 'pACBB_DmpR_BBa_J23112_BBa_B0032_L3S1P56', 'DW': 3.93}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A6', 'vol': 0.53, 'dil': 0, 'plate': 'pro'}, 'part3': {'well': 'B11', 'vol': 0.52, 'dil': 0, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p6_r23_t25', 'name': 'pACBB_DmpR_BBa_J23112_BBa_B0034_L3S1P56', 'DW': 5.51}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A12', 'vol': 1.9, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'A1', 'vol': 0.55, 'dil': 0, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p12_r1_t25', 'name': 'pACBB_DmpR_BBa_M13101_BBa_J61100_L3S1P56', 'DW': 4.11}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A12', 'vol': 1.9, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'B6', 'vol': 2.45, 'dil': 5, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p12_r18_t25', 'name': 'pACBB_DmpR_BBa_M13101_BBa_B0072_L3S1P56', 'DW': 2.21}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A12', 'vol': 1.9, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'B10', 'vol': 2.1, 'dil': 5, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p12_r22_t25', 'name': 'pACBB_DmpR_BBa_M13101_BBa_B0032_L3S1P56', 'DW': 2.56}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'A12', 'vol': 1.9, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'B11', 'vol': 0.52, 'dil': 0, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p12_r23_t25', 'name': 'pACBB_DmpR_BBa_M13101_BBa_B0034_L3S1P56', 'DW': 4.14}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'B6', 'vol': 1.75, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'A1', 'vol': 0.55, 'dil': 0, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p18_r1_t25', 'name': 'pACBB_DmpR_BBa_J23100_BBa_J61100_L3S1P56', 'DW': 4.26}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'B6', 'vol': 1.75, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'B6', 'vol': 2.45, 'dil': 5, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p18_r18_t25', 'name': 'pACBB_DmpR_BBa_J23100_BBa_B0072_L3S1P56', 'DW': 2.36}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'B6', 'vol': 1.75, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'B10', 'vol': 2.1, 'dil': 5, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p18_r22_t25', 'name': 'pACBB_DmpR_BBa_J23100_BBa_B0032_L3S1P56', 'DW': 2.71}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'B6', 'vol': 1.75, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'B11', 'vol': 0.52, 'dil': 0, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p18_r23_t25', 'name': 'pACBB_DmpR_BBa_J23100_BBa_B0034_L3S1P56', 'DW': 4.29}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'B9', 'vol': 2.4, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'A1', 'vol': 0.55, 'dil': 0, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p21_r1_t25', 'name': 'pACBB_DmpR_BBa_J23114_BBa_J61100_L3S1P56', 'DW': 3.61}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'B9', 'vol': 2.4, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'B6', 'vol': 2.45, 'dil': 5, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p21_r18_t25', 'name': 'pACBB_DmpR_BBa_J23114_BBa_B0072_L3S1P56', 'DW': 1.71}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'B9', 'vol': 2.4, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'B10', 'vol': 2.1, 'dil': 5, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p21_r22_t25', 'name': 'pACBB_DmpR_BBa_J23114_BBa_B0032_L3S1P56', 'DW': 2.06}}, {'part0': {'well': 'A1', 'vol': 5, 'dil': 0, 'plate': 'EXT'}, 'part1': {'well': 'A2', 'vol': 3.29, 'dil': 0, 'plate': 'EXT'}, 'part2': {'well': 'B9', 'vol': 2.4, 'dil': 5, 'plate': 'pro'}, 'part3': {'well': 'B11', 'vol': 0.52, 'dil': 0, 'plate': 'rbs'}, 'part4': {'well': 'C1', 'vol': 1.15, 'dil': 5, 'plate': 'ter'}, 'meta': {'No': 'pBB_Dm_p21_r23_t25', 'name': 'pACBB_DmpR_BBa_J23114_BBa_B0034_L3S1P56', 'DW': 3.64}}, ['ter', 'rbs', 'pro']]

def run(protocol: protocol_api.ProtocolContext):

    ## get global parameters
    well_data = meta_data[:-1]
    enz_vol = 4

    # Deck Setting
    ## Modules  
    module_thermocycler = protocol.load_module("thermocycler Module")

    ## Racks
    globals()['ter'] = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 1)
    globals()['rbs'] = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 2)
    globals()['pro'] = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 3)

    assemble_plate = module_thermocycler.load_labware("biorad_96_wellplate_200ul_pcr")
    EXT = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', 9)
    trash = protocol.loaded_labwares[12]["A1"]

    tiprack_20_1 = protocol.load_labware("opentrons_96_tiprack_20ul", 6)
    tiprack_20_2 = protocol.load_labware("opentrons_96_tiprack_20ul", 5)
    p20_sin = protocol.load_instrument("p20_single_gen2", "right", tip_racks=[tiprack_20_1, tiprack_20_2])

    ## Start Tiprack positions
    p20_sin.starting_tip = tiprack_20_1.well("A1")

    ## Reagents
    ### Every reagent should be in 1.5ml Bioneer screw tube.
    enz_mix = EXT['D1']
    DW = EXT['D5']


    # Protocol
    ## Functions
    def enzyme_transfer(pipette, volume, src, dest, delay_second=[0, 0],
                    top_delay=False, asp_rate=None, dis_rate=None,
                    mix_after=False, drop_tip = True):
        # top_delay = list, delay_second= list, mix_after = list

        if asp_rate:
            pipette.flow_rate.aspirate=asp_rate
        if dis_rate:
            pipette.flow_rate.dispense=dis_rate            

        if pipette._has_tip == False:
            pipette.pick_up_tip()

        pipette.aspirate(volume, src)
        protocol.delay(seconds=delay_second[0])
        if top_delay:
            pipette.move_to(src.top(z=-3))
            protocol.delay(seconds=top_delay[0])
        pipette.dispense(volume, dest)
        protocol.delay(seconds=delay_second[1])
        pipette.dispense(1, dest.top(z=-3))
        if type(mix_after) == list:
            try:
                pipette.flow_rate.aspirate=mix_after[2]
                pipette.flow_rate.dispense=mix_after[2]
            except:
                pass
            pipette.mix(mix_after[0], mix_after[1])
        if top_delay:
            pipette.move_to(src.top(z=-3))
            protocol.delay(seconds=top_delay[1])
        if drop_tip:    
            pipette.drop_tip()

    ## Enzyme Transfer
    ### get repeated_use value
    enz_dis_num = (19.5 / enz_vol)
    if enz_dis_num > 4:
        enz_dis_num = 4
    else:
        enz_dis_num = floor(enz_dis_num)

    ### Enzyme distribute
    src = enz_mix
    dest_wells = assemble_plate.wells()[:len(well_data)]

    n = 0
    while(1):
        if n + enz_dis_num > len(well_data):
            dest_tmp = dest_wells[n:]
        else:
            dest_tmp = dest_wells[n:n+enz_dis_num]

        if p20_sin._has_tip == False:
            p20_sin.pick_up_tip()
        p20_sin.flow_rate.aspirate = 3
        p20_sin.flow_rate.dispense = 3
        p20_sin.aspirate(enz_vol*enz_dis_num + 0.5, src)
        #delay
        protocol.delay(seconds=1)
        #topdelay
        p20_sin.move_to(src.top(z=-3))
        protocol.delay(seconds=1)

        for dest in dest_tmp:
            p20_sin.dispense(enz_vol, dest)
            #delay
            protocol.delay(seconds=1)

        p20_sin.drop_tip()
        n += enz_dis_num
        if n > enz_dis_num:
            break

    ### Part Transfer
    for i in range(well_data):
        data = well_data[i]
        dest = assemble_plate.wells()[i]

        for i2 in range(len(data)-1):
            tmp = data.get(f'part{i2}')
            vol = tmp['vol']
            # plate가 ext 포지션일 때는 name을 기준으로 well을 정해두기.
            src = eval(tmp['plate']).wells_by_name()[tmp['well']]
            enzyme_transfer(p20_sin, vol, src, dest,
                            asp_rate = 5, dis_rate =5, drop_tip=True)


    ## Thermocycling
    if module_thermocycler.lid_position == 'close':
        module_thermocycler.open_lid()

    module_thermocycler.set_temperature(37)
    module_thermocycler.set_temperature(37)
    module_thermocycler.set_temperature(37)

    protocol.pause("Protocol END \nIf you close this message, thermocycler open and deactivate.")
    protocol.disconnect()