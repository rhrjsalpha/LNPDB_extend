from kg_gen import KGGen

large_text = """
Lipid nanoparticles (LNPs) were formulated using ionizable lipid 304O13,
cholesterol, DSPC, and DMG-PEG2000. The molar ratio was 50:38.5:10:1.5.
These LNPs were tested in HeLa cells using Factor VII siRNA.
Dose-dependent knockdown of Factor VII was observed.
"""

kg = KGGen(
    model="ollama_chat/jinbora/deepseek-r1-Bllossom:8b",
    temperature=0.0
)

graph = kg.generate(
    input_data=large_text,
    chunk_size=2000,
    cluster=True,
    context="Lipid nanoparticle formulation and experiments"
)

print("Entities:", graph.entities)
print("Edges:", graph.edges)
print("Relations:", graph.relations)

KGGen.visualize(
    graph,
    output_path="kg_large.html",
    open_in_browser=True
)