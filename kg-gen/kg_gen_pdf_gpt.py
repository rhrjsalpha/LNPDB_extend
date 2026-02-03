from kg_gen import KGGen

large_text = """
The corresponding IL, 18:1 Δ9-cis phosphoethanolamine (DOPE),
cholesterol, and 14:0 PEG2000 phosphoethanolamine (C14-PEG2000)
were dissolved in ethanol at a molar ratio of 35:16:46.5:2.5, respectively,
to form the organic phase of the formulation. For the B10 formulation,
a molar ratio of 40:30:25:2.5 was used. The SM-102 LNP was
prepared using SM-102, 1,2-distearoyl-sn-glycero-3-phosphocholine
(DSPC), cholesterol, and 1,2-dimyristoyl-rac-glycero-3-methoxypolyethylene
glycol-2000 (DMG-PEG2000) at a molar ratio of
50:10:38.5:1.5, respectively. The ALC-0315 LNPwas formulated via ALC-
0315, DSPC, cholesterol, and ALC-0159 at a molar ratio of
46.3:9.4:42.7:1.6, respectively. The aqueous phase was prepared by
dissolving the corresponding mRNA in a 10mM citrate buffer at pH 3
(Teknova, Hollister, CA, USA). For each mRNA LNP, the weight ratio of
the mRNA to IL was 1:10, and the volume ratio of the organic phase to
aqueous phase was 1:3. Each phase was loaded in separate glass syringes
(Hamilton Company, Reno,NV) and connected to a Pump 33 DDS
syringe pump(HarvardApparatus,MA,USA) attached to amicrofluidic
device with a staggered herringbone micromixer design. The microfluidic
devices were fabricated in polydimethylsiloxane utilizing standard
soft lithographic procedures39. A two-step exposure process was
used to create the SU-8 master with positive channel features on a
silicon wafer, where each mixing channel is 4 cm in length. The syringes
were injected at a flow rate of 0.6mL/min and 1.8mL/min for the
organic phase and aqueous phase, respectively.
"""

kg = KGGen(
    model="openai/gpt-5.2",
    temperature=1.0,
    api_key=""
)

graph = kg.generate(
    input_data=large_text,
    chunk_size=2000,
    cluster=True,
    context="Lipid nanoparticle formulation and experiments"
)

print("Entities 1:", graph.entities)
print("Edges 1:", graph.edges)
print("Relations 1:", graph.relations)

#graph_clustered = kg.cluster(
#    graph,
#    context="Lipid nanoparticle formulation"
#)
#
#print("Entities 2:", graph.entities)
#print("Edges 2:", graph.edges)
#print("Relations 2:", graph.relations)

KGGen.visualize(
    graph,
    output_path="kg_large_gpt-5.2.html",
    open_in_browser=True
)