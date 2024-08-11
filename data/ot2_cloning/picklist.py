"""
Make picklist for liquid handler.
Using plateo [https://edinburgh-genome-foundry.github.io/Plateo/index.html].

Convert workflow table (stacked picklist) to unstacked piclist.
"""
import plateo
import pandas as pd

# example workflow table
workflow_df = pd.DataFrame({
    "name": ["frag1", "frag2"],
    "dna1": ["template1", "template2"],
    "dna2": ["primer1", "primer3"],
    "dna3": ["primer2", "primer4"],
    "dna4": ["[E]PCRmix", "[E]PCRmix"],
    "dna5": ["[E]dw", "[E]dw"],
})

