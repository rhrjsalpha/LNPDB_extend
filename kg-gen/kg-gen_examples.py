from kg_gen import KGGen

kg = KGGen(
    model="ollama_chat/jinbora/deepseek-r1-Bllossom:8b",   # 또는 ollama_chat/deepseek-r1:14b
    temperature=0.0
)

text = "Linda is Josh's mother. Ben is Josh's brother."

graph = kg.generate(input_data=text)

print(graph.entities)
print(graph.edges)
print(graph.relations)

KGGen.visualize(
    graph,
    output_path="kg.html",
    open_in_browser=True
)

