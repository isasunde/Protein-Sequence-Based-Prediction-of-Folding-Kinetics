"""Retrieve and save raw K-Prot API responses for UniProt IDs."""

import logging
import os
from pathlib import Path

from dotenv import load_dotenv

from Protein_Sequence_Based_Prediction_of_Folding_Kinetics import (
    get_endpoint,
    save_raw_response,
)

# Configuration
PROJECT_DIR = Path(__file__).resolve().parent
RAW_DIR = PROJECT_DIR / "data" / "raw"

UNIPROT_IDS = [
    "O54310", "O66474", "O75400", "P00171", "P00282",
    "P00523", "P00524", "P01053", "P02417", "P06876",
    "P02640", "P02966", "P03034", "P04049", "P06241",
    "P06654", "P06748", "P0ABE7", "P07107", "P07311",
    "P07751", "P0A9X9", "P0AA04", "P0AEZ7", "P0AFG6",
    "P0CG48", "P0DQD3", "P11182", "P11961", "P12931",
    "P13479", "P14621", "P15891", "P23370", "P23727",
    "P24821", "P31787", "P32081", "P38507", "P39476",
    "P40040", "P41016", "P45481", "P46937", "P48052",
    "P54274", "P56278", "P61024", "P62942", "P62993",
    "P68206", "Q13158", "Q23551", "Q51918", "Q53291",
    "Q56311", "Q8CGF7", "Q8ZUR6", "Q9NYB0", "Q9RA57",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)
logger = logging.getLogger(__name__)


def acquire_data():
    """Retrieve and save one raw API response per UniProt ID."""

    load_dotenv(PROJECT_DIR / ".env")

    api_url = os.getenv("API_URL")
    if not api_url:
        raise ValueError("API_URL is missing from the environment.")

    RAW_DIR.mkdir(parents=True, exist_ok=True)

    successful = []
    failed = []

    for uniprot_id in UNIPROT_IDS:
        filename = f"kpro_uniprot_{uniprot_id}.json"

        logger.info("Retrieving %s", uniprot_id)

        try:
            url, response = get_endpoint(
                id=uniprot_id,
                API_URL=api_url,
                endpoint="uniprot",
            )

            if response.status_code != 200:
                logger.warning(
                    "Request failed for %s: HTTP %s",
                    uniprot_id,
                    response.status_code,
                )
                failed.append(uniprot_id)
                continue

            response_data = response.json()

            save_raw_response(
                data=response_data,
                url=url,
                endpoint=f"UniProt ID: {uniprot_id}",
                filename=filename,
                RAW_DIR=RAW_DIR,
            )

            successful.append(uniprot_id)
            logger.info("Saved %s", filename)

        except Exception:
            # Log the exception and continue with the remaining IDs.
            logger.exception("Acquisition failed for %s", uniprot_id)
            failed.append(uniprot_id)

    logger.info(
        "Acquisition complete: %d successful, %d failed, %d total.",
        len(successful),
        len(failed),
        len(UNIPROT_IDS),
    )

    if failed:
        logger.warning("Failed UniProt IDs: %s", ", ".join(failed))

    return {
        "successful": successful,
        "failed": failed,
    }


if __name__ == "__main__":
    acquire_data()
