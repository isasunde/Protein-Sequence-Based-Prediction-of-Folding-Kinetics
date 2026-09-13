from sklearn.linear_model import LinearRegression, Ridge
from xgboost import XGBRegressor

def get_models():
    """
    Return the regression models and their preprocessing settings.

    Returns
    -------
    dict
        Dictionary containing model instances and whether feature
        scaling should be applied.
    """
    models = {
        "LinearRegression": {
            "model": LinearRegression(),
            "scale": True,
        },

        "Ridge": {
            "model": Ridge(
                alpha=10,
            ),
            "scale": True,
        },

        "XGBoost": {
            "model": XGBRegressor(
                n_estimators=250,
                max_depth=4,
                learning_rate=0.05,
                subsample=0.9,
                colsample_bytree=0.9,
                random_state=42,
                n_jobs=-1,
            ),
            "scale": False,
        },
    }

    return models
