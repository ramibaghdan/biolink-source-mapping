# Reference files the pipeline needs. Override on the command line:
#   make all HGNC=/path/to/hgnc_complete_set.txt
HGNC   ?= hgnc_complete_set.txt
MONDO  ?= data/mondo/mondo_nodes.tsv
LABELS ?= data/ontology_labels/labels.tsv
GENE   ?= NCBIGene:6416

.PHONY: all fetch build validate figure clean help

help:
	@echo "make fetch     download mondo_nodes.tsv and the OLS4 label cache"
	@echo "make build     run the pipeline, write nodes.csv and edges.csv"
	@echo "make validate  check biolink terms and naming coverage"
	@echo "make figure    redraw the demo subgraph"
	@echo "make all       build, validate, figure"
	@echo ""
	@echo "HGNC=$(HGNC)"

# Network. Run once. Labels need nodes.csv, so build first if it is missing.
fetch:
	python src/fetch.py mondo --out $(MONDO)
	python src/fetch.py labels --out $(LABELS)

build:
	python src/pipeline.py --hgnc $(HGNC) --mondo-nodes $(MONDO) --labels $(LABELS)

validate:
	python src/validate.py

figure:
	python src/visualize.py --gene $(GENE)

all: build validate figure

clean:
	rm -f output/nodes.csv output/edges.csv output/mondo_normalization.csv output/naming.csv
