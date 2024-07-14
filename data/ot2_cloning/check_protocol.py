def get_values(*names):
    import json
    _all_values = json.loads("""
                             export_JSON
                             """)
    return [_all_values[n] for n in names]

# Check needed tips
def check_tips():
    return 96
    pass

# Enzyme position
def enzyme_position(enzyme_list):
    # 24-well position
    enzyme_list = enzyme_list + ["DW", "CPcell", "LB"]
    well_position = ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4","C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"]
    return_dict = {}
    for enzyme, well in zip(enzyme_list, well_position):
        return_dict[enzyme] = well
    return return_dict

def deck_position(plates):
    position = [1,2,3,4,5,6,9]
    
    deck_dict = {}
    deck_dict["Enzyme_tube"] = position.pop(0)
    deck_dict["p20_tip"] = position.pop(0)
    deck_dict["p300_tip"] = position.pop(0)
    
    for key in plates.keys():
        if plates[key]["type"] == "Reaction":
            deck_dict[key] = 7
            continue
        
        deck_dict[key] = position.pop(0)
        
    return deck_dict