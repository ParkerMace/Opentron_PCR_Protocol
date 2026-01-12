from opentrons import protocol_api
import math

'''
The protocol for transformations to be tested with Sam. 
Main Points:
-Integrated with moclomatic
-Create agar wells?, calculate agar height based on volume
-Conserve pipette usage very conservatively due to space constraints
-Dilutions begin in one four square corner and move clockwise through the square
-add preheat for heatshock?
'''

metadata = {
    "protocol_name" : "Yeast_PCR",
    "author" : "Parker_Mace",
    'description': 'Automated PCR'
}

requirements = {
    "robotType" : "OT-2",
    "apiLevel" : "2.27"
}

def run(protocol: protocol_api.ProtocolContext):

    #---------
    # Labware
    # --------    

    #Modules
    hs_mod = protocol.load_module('heaterShakerModuleV1', '10')
    hs_adapter = hs_mod.load_adapter("opentrons_96_flat_bottom_adapter")
    recover = hs_adapter.load_labware("nest_96_wellplate_200ul_flat")

    temp_mod = protocol.load_module(module_name="temperature module gen2", location="1")
    temp_adapter = temp_mod.load_adapter("opentrons_96_well_aluminum_block")
    transform = temp_adapter.load_labware("nest_96_wellplate_100ul_pcr_full_skirt")

    #Resevoir and Assemblies
    resevoir = protocol.load_labware('opentrons_tough_4_reservoir_72ml', '2')
    assemblies = protocol.load_labware('opentrons_tough_4_reservoir_72ml', '5')

    #Pipette and tips
    p20_tiprack = protocol.load_labware('opentrons_96_tiprack_20ul', '4')
    p20_multi = protocol.load_instrument('p20_multi_gen2', 'left', tip_racks=[p20_tiprack])

    p300_tiprack = protocol.load_labware('opentrons_96_tiprack_300ul', '7')
    p300_multi = protocol.load_instrument('p300_multi_gen2', 'right', tip_racks=[p300_tiprack] )

    #Dilution Plate
    dil_plate = protocol.load_labware('appliedbiosystemsmicroamp_384_wellplate_40ul', '8')
    
    #Agar Plates
    agar_plate_1 = protocol.load_labware('axygen_96_wellplate_500ul', '3')
    agar_plate_2 = protocol.load_labware('axygen_96_wellplate_500ul', '6')
    agar_plate_3 = protocol.load_labware('axygen_96_wellplate_500ul', '9')

    #---------
    # Important Variables
    #--------  

    p20_multi.flow_rate.aspirate = 2.5
    p20_multi.flow_rate.dispense = 6

    p300_multi.flow_rate.aspirate = 40
    p300_multi.flow_rate.dispense = 70

    samples = 12
    sample_col = samples // 8
    assembly_wells = assemblies.wells()[0:sample_col:8]
    recovery_wells = recover.wells()[0:sample_col:8]
    trans_wells = transform.wells()[0:sample_col:8]

    plates = (agar_plate_1, agar_plate_2, agar_plate_3)

    dilution_wells = {}
    for i in range(sample_col):  # sample_col = number of 8-channel samples
        col_start = 2 * i

        dilution_wells[i] = (
            dil_plate.rows()[0][col_start],     # A?
            dil_plate.rows()[0][col_start + 1], # A?
            dil_plate.rows()[1][col_start],     # B?
            dil_plate.rows()[1][col_start + 1], # B?
        )
    
    agar_vol = 30
    agar_height = (agar_vol * 0.001) / (math.pi * math.sqrt(3.43))

    def create_plates(plate_vol):
        '''Create 96-well plates using tempered agar from resevoir.'''
        p300_multi.pick_up_tip()

        for i in range(sample_col):
            for plate in plates:
                dest = plate.columns()[i][0]
                p300_multi.aspirate(location=resevoir.wells()[1], volume = plate_vol)
                p300_multi.dispense(location = dest)
                p300_multi.blow_out(dest.top())
        p300_multi.return_tip()
        
    def distribute_media(dilution_vol, recovery_vol):
        """Distribute media to dilution and recovery wells."""
        p300_multi.pick_up_tip()
        for i in range(sample_col):
            for well in dilution_wells[i]:
                p300_multi.aspirate(location=resevoir.wells()[0], volume = dilution_vol)
                p300_multi.dispense(location = well)
                p300_multi.blow_out(well.top())
        for well in recovery_wells:
            p300_multi.aspirate(location=resevoir.wells()[0], volume = recovery_vol)
            p300_multi.dispense(location = well)
            p300_multi.blow_out(well.top())
        p300_multi.return_tip()

    def transformation(assembly_vol):
        '''
        Uses heater/shaker to transform cells
        
        :param assembly_vol: Volume to be transformed per well.
        '''
        for well in assembly_wells:
            dest = transform.well(well.well_name)
            p20_multi.pick_up_tip()
            p20_multi.aspirate(location= well, volume = assembly_vol)
            p20_multi.dispense(location = dest)
            p20_multi.blow_out(well.top())
            p20_multi.return_tip()

    def recovery(recover_vol):
        '''
        Uses temp module to recover transformed cells
        
        :param trans_vol: Volume to be recovered per well.
        '''
        for well in trans_wells:
            dest = recover.well(well.well_name)
            p20_multi.pick_up_tip()
            p20_multi.aspirate(location= well, volume = recover_vol)
            p20_multi.dispense(location = dest)
            p20_multi.blow_out(well.top())
            p20_multi.return_tip()

    def dilutions(dil1,dil2,dil3,dil4):
        '''Dilution of transformed samples.'''

        for i in range(sample_col):
            sample = recover.columns()[i][0]
            '''
            Dilution order:
            [1|3]
            [2|4]
            '''
            p20_multi.pick_up_tip()
            p20_multi.aspirate(location = sample, volume = dil1)
            p20_multi.dispense(location = dilution_wells[i][0])
            p20_multi.mix(repetitions = 5, volume=3)
            p20_multi.aspirate(location = dilution_wells[i][0], volume = dil2)
            p20_multi.dispense(location = dilution_wells[i][1])
            p20_multi.mix(repetitions = 5, volume=3)
            p20_multi.aspirate(location = dilution_wells[i][1], volume = dil3)
            p20_multi.dispense(location = dilution_wells[i][2])
            p20_multi.mix(repetitions = 5, volume=3)
            p20_multi.aspirate(location = dilution_wells[i][2], volume = dil4)
            p20_multi.dispense(location = dilution_wells[i][3])
            p20_multi.mix(repetitions = 5, volume=3)
            p20_multi.return_tip()

    def plating():
        '''
        Plating of dilutions on 96-well agar plates.
        '''
        p20_multi.well_bottom_clearance.dispense = agar_height + 0.3

        for i in range(sample_col):
                p20_multi.pick_up_tip()
                p20_multi.aspirate(location = dilution_wells[i][1], volume = 3)
                p20_multi.dispense(location = agar_plate_1.wells()[i])
                p20_multi.aspirate(location = dilution_wells[i][2], volume = 3)
                p20_multi.dispense(location = agar_plate_2.wells()[i])
                p20_multi.aspirate(location = dilution_wells[i][3], volume = 3)
                p20_multi.dispense(location = agar_plate_3.wells()[i])
                p20_multi.return_tip()

    def main():
        hs_mod.close_labware_latch()

        #create agar plates
        create_plates(30)

        #on ice
        on_ice = temp_mod.start_set_temperature(celsius=4)

        #preheat heater/shaker
        hs_mod.set_target_temperature(celsius=37)

        #Distribute media into recovery and dilution wells
        distribute_media(60, 90)

        # transformation profile
        protocol.wait_for_tasks([on_ice])
        transformation(9)
        temp_mod.set_temperature(celsius=40)
        protocol.delay(seconds = 30)
        temp_mod.set_temperature(celsius=4)
        protocol.delay(minutes = 10)
        temp_mod.deactivate()

        #recovery profile
        recovery(20)
        hs_mod.set_and_wait_for_temperature(37)
        hs_mod.set_and_wait_for_shake_speed(300)
        protocol.delay(minutes = 60)
        hs_mod.deactivate_heater()
        hs_mod.deactivate_shaker()

        #create dilutions
        dilutions(3,3,3,3)
        hs_mod.deactivate_heater()

        #plate dilutions
        plating()
        

    #run protocol
    main()

        
