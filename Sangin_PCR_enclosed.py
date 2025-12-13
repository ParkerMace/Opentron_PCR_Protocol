from opentrons import protocol_api
import math
from itertools import islice

#----------------------------------------
#PCR for Double and Triple Gene Knockouts - Auxin Synthesis Pathway
#----------------------------------------

metadata = {
    "protocol_name" : "Yeast_PCR",
    "author" : "Parker_Mace",
    'description': 'Automated PCR'
}

requirements = {
    "robotType" : "OT-2",
    "apiLevel" : "2.16"
}

def run(protocol: protocol_api.ProtocolContext):

    # ----------------------
    # USER PARAMETERS (edit volumes as needed)
    # ----------------------
    replicates = 1
    water_per_rxn = 21  # µL
    onetaq_per_rxn = 25  # µL
    overage = 1.1
    vol_master_mix = water_per_rxn + onetaq_per_rxn # uL per well
    vol_primer = 2 # uL per primer 
    vol_dna = 2 # uL per dna sample
    vol_reaction = vol_primer + vol_dna + vol_master_mix

    

    # ----------------------
    # LABWARE (adjust labware and deck positions to your physical setup)
    # ----------------------

    # Load the Heater-Shaker with your PCR plate
    tc_mod = protocol.load_module(module_name="thermocycler")
    pcr_plate = tc_mod.load_labware(name='opentrons_96_wellplate_200ul_pcr_full_skirt')

    p300_tiprack = protocol.load_labware('opentrons_96_tiprack_300ul', '6')
    p300 = protocol.load_instrument('p300_single_gen2', 'left', tip_racks=[p300_tiprack])

    p20_tiprack = protocol.load_labware('opentrons_96_tiprack_20ul', '3')
    p20 = protocol.load_instrument('p20_single_gen2', 'right', tip_racks=[p20_tiprack])

    master_mix_tuberack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', '2')
    master_mix_tube = master_mix_tuberack.wells_by_name()['A1']
    water = master_mix_tuberack.wells_by_name()['A2']
    onetaq = master_mix_tuberack.wells_by_name()['A3']

    primer_rack = protocol.load_labware('opentrons_96_aluminumblock_generic_pcr_strip_200ul', '5')

    dna_plate = protocol.load_labware('opentrons_96_aluminumblock_generic_pcr_strip_200ul', '4')

    '''
    Use PCR_OT2_CSV_to_dict to convert csv to dictionary formats and paste below.
    Or input manually in correct format.
    '''
    
    # Copy/Paste dictionaries from generated text file
    primer_map = {"ARO8": primer_rack.wells_by_name()["A1"],
    "NIT1": primer_rack.wells_by_name()["B1"],
    "AMD": primer_rack.wells_by_name()["C1"],}
    sample_genes = {"M29": ["ARO8", "NIT1","AMD"],}
    

    '''Other important dictionaries and lists.'''
    dna_sources = {}
    for i, sample in enumerate(sample_genes.keys()):
        dna_sources[sample] = dna_plate.wells()[i]

    protocol.comment("DNA sample sources mapped:")
    for s, w in dna_sources.items():
        protocol.comment(f"  {s} -> {w.display_name}")

    # --- list of all reactions (sample,gene,replicate) This will result in samples and replicates grouped together. ---
    reaction_list = []
    for sample, genes in sample_genes.items():
        if not genes:
            continue
        for r in range(replicates):
            for gene in genes:
                reaction_list.append((sample, gene, r))

    reaction_assignments = {}
    for dest, (sample, gene, replicate) in zip(pcr_plate.wells(), reaction_list):
        primer_well = primer_map[gene]
        reaction_assignments[dest] = (sample, gene, replicate, primer_well)

    total_reactions = len(reaction_list)
    if total_reactions == 0:
        raise RuntimeError("No reactions to plate (no sample–gene combos found).")

    plates_needed = math.ceil(total_reactions / 96)
    protocol.comment(f"Total reactions: {total_reactions}. Plates required: {plates_needed}.")

    # ----------------------
    # Important functions
    # ----------------------
    def create_master_mix():
        """Create master mix with overage for total reactions."""
        water_vol = water_per_rxn * total_reactions * overage
        onetaq_vol = onetaq_per_rxn * total_reactions * overage

        protocol.comment(f"Creating master mix: {water_vol:.1f} µL water + {onetaq_vol:.1f} µL OneTaq")

        max_vol = p300.max_volume
        
        # Transfer water
        p300.pick_up_tip()
        remaining = water_vol
        while remaining > 0:
            transfer_vol = min(remaining, max_vol)
            p300.aspirate(transfer_vol, water)
            p300.dispense(transfer_vol, master_mix_tube)
            p300.blow_out(master_mix_tube.top())
            remaining -= transfer_vol
        p300.drop_tip()

        # Transfer OneTaq
        p300.pick_up_tip()
        remaining = onetaq_vol
        while remaining > 0:
            transfer_vol = min(remaining, max_vol)
            p300.aspirate(transfer_vol, onetaq)
            p300.dispense(transfer_vol, master_mix_tube)
            p300.blow_out(master_mix_tube.top())
            remaining -= transfer_vol
        
        # Mix master mix
        p300.mix(5, 200, master_mix_tube)
        p300.blow_out(master_mix_tube.top())
        p300.drop_tip()

        protocol.comment("Master mix prepared and mixed.")

    def distribute_master_mix(dest_wells):
        """Distribute master mix to wells, one at a time."""
        p300.pick_up_tip()
        for well in dest_wells:
            p300.aspirate(vol_master_mix, master_mix_tube.bottom(2))
            p300.dispense(vol_master_mix, well.bottom(2))
            p300.blow_out(well.top())
        p300.drop_tip()


    def add_primers(dest_wells, reaction_assignments):
        """Add primers to destination wells."""
        for dest in dest_wells:
            sample, gene, replicate, primer_well = reaction_assignments[dest]
            p20.pick_up_tip()
            p20.aspirate(vol_primer, primer_well, rate=0.5)
            p20.dispense(vol_primer, dest, rate=0.5)
            p20.mix(3, 10, dest)
            p20.blow_out(dest.top())
            p20.touch_tip()
            p20.drop_tip()

    def add_dna(dest_wells, reaction_assignments, dna_sources):
        """Add DNA samples to destination wells."""
        for dest in dest_wells:
            sample, gene, replicate, primer_well = reaction_assignments[dest]
            dna_source = dna_sources[sample]
            p20.pick_up_tip()
            p20.aspirate(vol_dna, dna_source, rate=0.5)
            p20.dispense(vol_dna, dest, rate=0.5)
            p20.mix(5, 15, dest)
            p20.blow_out(dest.top())
            p20.touch_tip()
            p20.drop_tip()

    def run_pcr(
            denature_temp: float,
            denature_time: int,
            anneal_temp: float,
            anneal_time: int,
            extend_temp: float,
            extend_time: int,
            elongation_step_temp: float,
            elongation_step_time: int,
            cycles: int,
            final_hold: float = 4.0
        ):
            """
            Run a PCR cycle on the Heater-Shaker module.

            Args:
                denature_temp (°C), denature_time (sec)
                anneal_temp (°C), anneal_time (sec)
                extend_temp (°C), extend_time (sec)
                cycles (int): number of PCR cycles
                final_hold (°C): temperature to hold at the end
            """

            protocol.comment(f"Starting PCR program: {cycles} cycles")

            tc_mod.close_lid()
            tc_mod.set_lid_temperature(temperature = 105)

            for cycle in range(1, cycles + 1):
                protocol.comment(f"Cycle {cycle} / {cycles}")

                # Denaturation
                tc_mod.set_block_temperature(
                temperature=denature_temp,
                hold_time_seconds=denature_time,
                block_max_volume=vol_reaction)

                # Annealing
                tc_mod.set_block_temperature(
                temperature=anneal_temp,
                hold_time_seconds=anneal_time,
                block_max_volume=vol_reaction)

                # Extension
                tc_mod.set_block_temperature(
                temperature=extend_temp,
                hold_time_seconds=extend_time,
                block_max_volume=vol_reaction)

            # Final Elongation
            tc_mod.set_block_temperature(
                temperature=elongation_step_temp,
                hold_time_seconds=elongation_step_time,
                block_max_volume=vol_reaction)
            # Final hold
            tc_mod.deactivate_lid()
            tc_mod.set_block_temperature(temperature=final_hold, block_max_volume=vol_reaction)
            protocol.set_rail_lights(True)
            protocol.comment(f"PCR complete. Holding at {final_hold} °C.")
    
    # helper to chunk the reactions for each plate
    def chunked_iterable(iterable, size):
        it = iter(iterable)
        while True:
            chunk = list(islice(it, size))
            if not chunk:
                break
            yield chunk

    protocol.pause("Ensure reagents are loaded: water in A2, OneTaq in A3, empty 1.5mL tube in A1.")
    
    # Create master mix once for all reactions
    create_master_mix()

    # For each plate, ask user to place an empty PCR plate in the configured plate slot (same slot reused)
    plate_number = 1
    reaction_iter = chunked_iterable(reaction_list, 96)
    for plate_chunk in reaction_iter:
        protocol.comment(f"=== Preparing plate {plate_number} of {plates_needed} (contains {len(plate_chunk)} reactions) ===")

        protocol.pause("Ensure correct amount of tubes and reagents are placed in the modules.")

        # build the target well list for this plate (wells in order A1,A2...H12)
        target_wells = [pcr_plate.wells()[i] for i in range(len(plate_chunk))]

        # Open thermocycler
        tc_mod.open_lid()

        # 1) Distribute master mix
        protocol.comment("Adding master mix...")
        distribute_master_mix(target_wells)

        # 2) Add primers
        protocol.comment("Adding primers...")
        add_primers(target_wells, reaction_assignments)

        # 3) Add DNA
        protocol.comment("Adding DNA samples...")
        add_dna(target_wells, reaction_assignments, dna_sources)

        # 4) Run PCR
        protocol.pause("Cap PCR tubes.")
        run_pcr(
            denature_temp=94,
            denature_time=90,
            anneal_temp=56,
            anneal_time=45,
            extend_temp=59,
            extend_time=60,
            elongation_step_temp=68,
            elongation_step_time=300,
            cycles=30,
            final_hold=4
            )
        
        # 5) Completion
        tc_mod.open_lid()
        protocol.pause(f"Plate {plate_number} complete. Remove samples and press Resume for next plate.")

        # 6) Move on to next batch.
        plate_number += 1
