PY=python src

.PHONY: reproduce data preprocess features classify biomarkers pathways extensions sensitivity figure test notebook clean

## full pipeline, steps 0-7 (run `make data` first)
reproduce: preprocess features classify biomarkers pathways extensions sensitivity figure
	@echo "Pipeline complete. See results/ and docs/."

data:
	bash download_data.sh

preprocess:
	$(PY)/00_preprocess.py

features:
	$(PY)/01_feature_selection.py

classify:
	$(PY)/02_classification.py

biomarkers:
	$(PY)/03_biomarkers_roc.py

pathways:
	$(PY)/04_pathways.py

extensions:
	$(PY)/05a_three_class.py
	$(PY)/05b_shap_ablation.py
	$(PY)/05c_survival.py
	$(PY)/05d_novel.py

sensitivity:
	$(PY)/06_sensitivity.py

figure:
	$(PY)/07_summary_figure.py

## smoke test on synthetic data (no downloads, about a minute)
test:
	pytest -q tests

## re-execute the walkthrough notebook against results/
notebook:
	jupyter nbconvert --to notebook --execute --inplace notebooks/walkthrough.ipynb

# stage caches in data/processed make reruns resumable; clean forces a full recompute
clean:
	rm -rf data/processed/* results/figures/* results/tables/*
