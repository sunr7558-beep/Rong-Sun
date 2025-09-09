- Situation: Existing cardiac MRI models had poor accuracy distinguishing normal vs. abnormal scans.
- Task: Reformulate into a robust binary classifier with reproducibility and interpretability.
- Action: Built preprocessing pipeline (slice extraction, normalization), implemented CNNs in Python, applied class re-weighting & augmentation, added Grad-CAM saliency maps.
- Result: Accuracy improved from 74% → 86%, ROC-AUC ≈ 0.865, F1 ≈ 0.914; saliency maps highlighted clinically relevant regions, improving clinician trust.


**Data**: The dataset can be downloaded from [Kaggle link](https://www.kaggle.com/datasets/samdazel/automated-cardiac-diagnosis-challenge-miccai17?resource=download).
