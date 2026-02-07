from kg_gen import KGGen

large_text = """
1. An ionizable lipid of Formula A:
Figure US20170210697A1-20170727-C00303
or a salt thereof.
2. The ionizable lipid of claim 1, wherein the salt is a pharmaceutically acceptable salt.
3. A nanoparticle composition comprising a lipid component comprising an ionizable lipid of claim 1.
4. The nanoparticle composition of claim 3, wherein the lipid component further comprises a phospholipid.
5. The nanoparticle composition of claim 4, wherein the phospholipid is selected from the group consisting of 1,2-dilinoleoyl-sn-glycero-3-phosphocholine (DLPC),
1,2-dimyristoyl-sn-glycero-phosphocholine (DMPC), 1,2-dioleoyl-sn-glycero-3-phosphocholine (DOPC), 1,2-dipalmitoyl-sn-glycero-3-phosphocholine (DPPC),
1,2-distearoyl-sn-glycero-3-phosphocholine (DSPC), 1,2-diundecanoyl-sn-glycero-phosphocholine (DUPC), 1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine (POPC),
1,2-di-O-octadecenyl-sn-glycero-3-phosphocholine (18:0 Diether PC),
1-oleoyl-2-cholesterylhemisuccinoyl-sn-glycero-3-phosphocholine (OChemsPC),
1-hexadecyl-sn-glycero-3-phosphocholine (C 16 Lyso PC),
1,2-dilinolenoyl-sn-glycero-3-phosphocholine, 1,2-diarachidonoyl-sn-glycero-3-phosphocholine,
1,2-didocosahexaenoyl-sn-glycero-3-phosphocholine, 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine (DOPE), 1,2-diphytanoyl-sn-glycero-3-phosphoethanolamine (ME 16.0 PE),
1,2-distearoyl-sn-glycero-3-phosphoethanolamine,
1,2-dilinoleoyl-sn-glycero-3-phosphoethanolamine,
1,2-dilinolenoyl-sn-glycero-3-phosphoethanolamine,
1,2-diarachidonoyl-sn-glycero-3-phosphoethanolamine,
1,2-didocosahexaenoyl-sn-glycero-3-phosphoethanolamine,
1,2-dioleoyl-sn-glycero-3-phospho-rac-(1-glycerol) sodium salt (DOPG), sphingomyelin, and mixtures thereof.
6. The nanoparticle composition of claim 5, wherein the phospholipid is DSPC or DOPE.
7. The nanoparticle composition of claim 4, wherein the lipid component further comprises a structural lipid.
8. The nanoparticle composition of claim 7, wherein the structural lipid is selected from the group consisting of cholesterol, fecosterol, sitosterol, ergosterol, campesterol, stigmasterol, brassicasterol, tomatidine, ursolic acid, alpha-tocopherol, and mixtures thereof.
9. The nanoparticle composition of claim 8, wherein the structural lipid is cholesterol.
10. The nanoparticle composition of claim 7, wherein the lipid component further comprises a PEG lipid.
11. The nanoparticle composition of claim 10, wherein the PEG lipid is selected from the group consisting of a PEG-modified phosphatidylethanolamine, a PEG-modified phosphatidic acid, a PEG-modified ceramide, a PEG-modified dialkylamine, a PEG-modified diacylglycerol, a PEG-modified dialkylglycerol, and mixtures thereof.
12. The nanoparticle composition of claim 10, wherein the lipid component comprises about 30 mol % to about 60 mol % said ionizable lipid, about 0 mol % to about 30 mol % phospholipid, about 18.5 mol % to about 48.5 mol % structural lipid, and about 0 mol % to about 10 mol % PEG lipid.
13. The nanoparticle composition of claim 10, wherein the lipid component comprises about 35 mol % to about 55 mol % said ionizable lipid, about 5 mol % to about 25 mol % phospholipid, about 30 mol % to about 40 mol % structural lipid, and about 0 mol % to about 10 mol % PEG lipid.
14. The nanoparticle composition of claim 10, wherein the lipid component comprises about 50 mol % said ionizable lipid, about 10 mol % phospholipid, about 38.5 mol % structural lipid, and about 1.5 mol % PEG lipid.
15. The nanoparticle composition of claim 10, further comprising a therapeutic and/or prophylactic agent.
16. The nanoparticle composition of claim 12, further comprising a therapeutic and/or prophylactic agent.
17. The nanoparticle composition of claim 13, further comprising a therapeutic and/or prophylactic agent.
18. The nanoparticle composition of claim 15, wherein the therapeutic and/or prophylactic agent is a nucleic acid.
19. The nanoparticle composition of claim 15, wherein the therapeutic and/or prophylactic agent is a ribonucleic acid (RNA).
20. The nanoparticle composition of claim 19, wherein the RNA is selected from the group consisting of a small interfering RNA (siRNA), an asymmetrical interfering RNA (aiRNA), a microRNA (miRNA), a Dicer-substrate RNA (dsRNA), a small hairpin RNA (shRNA), a messenger RNA (mRNA), and mixtures thereof.
21. The nanoparticle composition of claim 19, wherein the RNA is an mRNA.
22. The nanoparticle composition of claim 21, wherein the mRNA includes one or more of a stem loop, a chain terminating nucleoside, a polyA sequence, a polyadenylation signal, and/or a 5′ cap structure.
23. The nanoparticle composition of claim 21, wherein the encapsulation efficiency of the therapeutic and/or prophylactic agent is at least 80% or at least 90%.
24. The nanoparticle composition of claim 21, wherein the wt/wt ratio of the lipid component to the mRNA is from about 10:1 to about 60:1.
25. The nanoparticle composition of claim 21, wherein the wt/wt ratio of the lipid component to the mRNA is about 20:1.
26. The nanoparticle composition of claim 21, wherein the N:P ratio is from about 5:1 to about 8:1.
27. A pharmaceutical composition comprising the nanoparticle composition of claim 21 and a pharmaceutically acceptable carrier.
28. A method of delivering a therapeutic and/or prophylactic agent to a mammalian cell, the method comprising administering to a subject the nanoparticle composition of claim 15, said administering comprising contacting the cell with the nanoparticle composition, whereby the therapeutic and/or prophylactic agent is delivered to the cell.
29. A method of producing a polypeptide of interest in a mammalian cell, the method comprising contacting the cell with the nanoparticle composition of claim 21, wherein the mRNA encodes the polypeptide of interest, whereby the mRNA is capable of being translated in the cell to produce the polypeptide of interest.
30. A method of synthesizing an ionizable lipid of claim 1, comprising reacting heptadecan-9-yl 8-((2-hydroxyethyl)amino)octanoate with nonyl-8-bromooctanoate under a suitable condition to provide the ionizable lipid of claim 1.
"""
print("Using model...")
kg = KGGen(
    model="openai/gpt-5.2",
    temperature=1.0,
    api_key=""
)

claim_context = """
You are extracting a PATENT CLAIM-LEVEL knowledge graph.

CRITICAL RULES:
1. Use ONLY information explicitly stated in the claims.
2. Ignore description, embodiments, examples, and Markush scaffolds.
3. Do NOT infer chemical structures beyond what is explicitly claimed.
4. Treat each independent claim as a primary entity.
5. Dependent claims only add constraints to their parent claim.
6. All entities MUST be traceable to one or more claim numbers.

ENTITY TYPES (use only these):
- ClaimedCompound
- ClaimedComposition
- ClaimedMethod
- Component
- Ratio
- Property
- Process
- Use

When molar ratios are given, decompose them into component-wise quantitative relations.
"""


print("Generating knowledge graph...")
graph = kg.generate(
    input_data=large_text,
    chunk_size=2000,
    cluster=True,
    context=claim_context
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

print("Visualizing knowledge graph...")
KGGen.visualize(
    graph,
    output_path="kg_gpt-5.2_US20170210697A1.html",
    open_in_browser=True
)