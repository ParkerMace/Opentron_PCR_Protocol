from opentrons import protocol_api
import math
from itertools import islice

#----------------------------------------
#PCR for Double and Triple Gene Knockouts - Auxin Synthesis Pathway
#----------------------------------------

metadata = {
    "protocol_name" : "April_Auxin_PCR",
    "author" : "Parker_Mace",
    'description': 'Automated PCR setup with master mix, primers, and DNA samples'
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
    vol_master_mix = 23.0 # uL per well
    vol_primer = 1 # uL per primer (6 primers => 3 uL total primers)
    vol_dna = 1 # uL per dna sample

    # ----------------------
    # LABWARE (adjust labware and deck positions to your physical setup)
    # ----------------------

    # Load the Heater-Shaker with your PCR plate
    tc_mod = protocol.load_module(module_name="thermocyclerModuleV1")
    pcr_plate = tc_mod.load_labware(name='opentrons_96_aluminumblock_generic_pcr_strip_200ul')

    p300_tiprack = protocol.load_labware('opentrons_96_tiprack_300ul', '6')
    p300 = protocol.load_instrument('p300_single_gen2', 'left', tip_racks=[p300_tiprack])

    p20_tiprack = protocol.load_labware('opentrons_96_tiprack_20ul', '3')
    p20 = protocol.load_instrument('p20_single_gen2', 'right', tip_racks=[p20_tiprack])

    master_mix_tuberack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', '2')
    master_mix_tube = master_mix_tuberack.wells_by_name()['A1']

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
        primer_well = primer_map[gene]   # or however your primers are stored
        reaction_assignments[dest] = (sample, gene, replicate, primer_well)

    total_reactions = len(reaction_list)
    if total_reactions == 0:
        raise RuntimeError("No reactions to plate (no sample–gene combos found).")

    plates_needed = math.ceil(total_reactions / 96)
    protocol.comment(f"Total reactions: {total_reactions}. Plates required: {plates_needed}.")

    # ----------------------
    # Important functions
    # ----------------------
    
    def distribute_master_mix(dest_wells):
        """
        Use one tip to distribute master mix to multiple wells,
        refilling from source tube as needed.
        """
        max_vol = p300.max_volume  # typically 300 uL
        per_well = vol_master_mix
        wells_remaining = list(dest_wells)

        p300.pick_up_tip()
        while wells_remaining:
            # How many wells can we serve in one aspiration?
            num_wells_this_round = int(max_vol // per_well)
            wells_chunk = wells_remaining[:num_wells_this_round]
            wells_remaining = wells_remaining[num_wells_this_round:]

            # Total volume to aspirate
            total_asp = per_well * len(wells_chunk)
            p300.aspirate(total_asp, master_mix_tube)

            # Dispense per well
            for w in wells_chunk:
                p300.dispense(per_well, w)

            # Small blowout back into source to clear tip
            p300.blow_out(master_mix_tube)

        p300.drop_tip()

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

            for cycle in range(1, cycles + 1):
                protocol.comment(f"Cycle {cycle} / {cycles}")

                # Denaturation
                tc_mod.set_block_temperature(
                temperature=denature_temp,
                hold_time_seconds=denature_time)

                # Annealing
                tc_mod.set_block_temperature(
                temperature=anneal_temp,
                hold_time_seconds=anneal_time)

                # Extension
                tc_mod.set_block_temperature(
                temperature=extend_temp,
                hold_time_seconds=extend_time)

            # Final Elongation
            tc_mod.set_block_temperature(
                temperature=elongation_step_temp,
                hold_time_seconds=elongation_step_time)
            # Final hold
            tc_mod.set_block_temperature(final_hold)
            protocol.comment(f"PCR complete. Holding at {final_hold} °C.")
    
    # helper to chunk the reactions for each plate
    def chunked_iterable(iterable, size):
        it = iter(iterable)
        while True:
            chunk = list(islice(it, size))
            if not chunk:
                break
            yield chunk

    # For each plate, ask user to place an empty PCR plate in the configured plate slot (same slot reused)
    plate_number = 1
    reaction_iter = chunked_iterable(reaction_list, 96)
    for plate_chunk in reaction_iter:
        protocol.comment(f"=== Preparing plate {plate_number} of {plates_needed} (contains {len(plate_chunk)} reactions) ===")

        # load/assume the plate is in the same slot (reuse plate slot). If you want multiple plates preloaded,
        # you can instead load multiple plate labware objects at the top and index them here.
        #pcr_plate = protocol.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', '5')

        # build the target well list for this plate (wells in order A1..H12)
        target_wells = [pcr_plate.wells()[i] for i in range(len(plate_chunk))]

        tc_mod.open_lid()

        # 1) Add master mix to all target wells using your multi-dispense helper
        distribute_master_mix(target_wells)

        # 2) Add primer mixes for this plate
        for dest, (sample, gene, replicate, primer_well) in reaction_assignments.items():
            if dest in target_wells:
                p20.transfer(vol_primer, primer_well, dest, new_tip='never')
        p20.drop_tip()

        # 3) Add DNA for this plate
        for dest, (sample, gene, replicate, primer_well) in reaction_assignments.items():
            if dest in target_wells:
                dna_source = dna_sources[sample]
                p20.transfer(vol_dna, dna_source, dest, new_tip='always')

        # 4) print a human-readable mapping for this plate (include primer mix well for clarity)
        protocol.comment(f"Plate {plate_number} mapping (first column = well on plate, sample, gene, replicate):")
        for well_idx, (sample, gene, rep_idx) in enumerate(plate_chunk):
            dest = target_wells[well_idx]
            primer_well = primer_map.get(gene)
            primer_loc = primer_well.display_name if primer_well is not None else "UNKNOWN"
            protocol.comment(f"{dest.display_name}: {sample} | {gene} | replicate {rep_idx+1} | primer mix: {primer_loc}")

        # 5) Run PCR
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

        plate_number += 1
