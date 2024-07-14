def flow_rate(pipette, **kwargs):
    assert (item in ['aspirate', 'dispense', 'blow_out'] for item in kwargs.keys()), "Error Keywords in Flow Rate."
    for i in kwargs.keys:
        pipette.flow_rate[i] = kwargs[i]

def enzyme_transfer(pipette, volume, src, dest, delay_second=[0, 0],
                top_delay=False, asp_rate=None, dis_rate=None,
                mix_after: list = [], drop_tip = True):
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

def get_well(well):
    #check가 되도록해야할듯.
    if type(well) == list:
        return (well[0])
    else:
        return (well)

def well_mix(pipette, dest, mix_params, flow_rate):
    if pipette._has_tip != True:
        pipette.pick_up_tip()
    pipette.flow_rate.aspirate = flow_rate
    pipette.flow_rate.dispense = flow_rate
    pipette.move_to(dest.bottom())
    pipette.mix(mix_params[0], mix_params[1])
    pipette.blow_out()
    pipette.move_to(dest.top(z=-3))
    protocol.delay(seconds=2)
    pipette.drop_tip()

def serial_distribute(pipette, volume, src, dest, asp_rate, dis_rate, max_repeat, delay_seconds=[0,0]):
    # dest should be well list
    def get_repeat_num():
        dis_num = (pipette.max_volume - 0.5) / volume
        if dis_num > max_repeat:
            return (max_repeat)
        else:
            return floor(dis_num)

    dis_num = get_repeat_num()
    n = 0

    if n + dis_num > len(dest):
        dest_tmp = dest[n:]
    else:
        dest_tmp = dest[n:n+dis_num]
    
    if pipette._has_tip == False:
        pipette.pick_up_tip()
    
    if asp_rate:
        pipette.flow_rate.aspirate = asp_rate
    if dis_rate:
        pipette.flow_rate.dispense = dis_rate
    
    pipette.aspirate(volume*dis_num+0.5, src)
    protocol.delay(seconds=delay_seconds[0])
    if top_delay:
        pass

    for dest_well in dest_tmp:
        pipette.dispense(volume, dest_well)
        #delay
        protocol.delay(seconds=delay_seconds[1])
    
    if drop_tip:
        pipette.drop_tip()

