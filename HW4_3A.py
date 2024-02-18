{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "a2518827-54a7-4fe8-8ba4-370a1b7e5319",
   "metadata": {},
   "source": [
    "A. Restriction enzymes are proteins that cleave DNA at specific sequence motifs. Use the file posted under week 4 materials, pGL3.fa(a linearized plasmid sequence file) "
   ]
  },
  {
   "cell_type": "markdown",
   "id": "0defe962-5576-47de-81b7-47470e557bf7",
   "metadata": {},
   "source": [
    "1. Write python code to read the plasmid sequence and returns the number and locations of matches to each the following restriction enzymes:"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "1641eb88-18c7-4cb8-ace4-0c1e8114bd45",
   "metadata": {},
   "source": [
    "ApoI (motif 5' RAATTY 3'; note that R = A or G, Y = T or C), BsaI (motif 5' GGTCTC 3'), and DpnI (GATC)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "ee5905e6-4c7a-4e15-8de6-c0544dfcd15b",
   "metadata": {},
   "source": [
    "Your script should search for both positive (sense) and negative (antisense) strand matches (NOTE: not all sites are palindromes).\n",
    "Finally, include code to digest with all three enzymes simultaneously"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "4ff3c29c-5e19-40ca-8532-1411f976d15b",
   "metadata": {},
   "source": [
    "This criterion is linked to a Learning OutcomeRestriction digest script - reads plasmid sequence file\n",
    "Code should open a filehandle and read lines. This can use pyfaidx or biopython seqIO, or can be manually written"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "6b3d21ce-04a8-4cd3-8e5d-d1a7d65d4218",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Results for ApoI:\n",
      "Site: 562\n",
      "Site: 1918\n",
      "Site: 4556\n",
      "Site: 4567\n",
      "Results for BsaI:\n",
      "Site: 3285\n",
      "Results for DpnI:\n",
      "Site: 109\n",
      "Site: 116\n",
      "Site: 741\n",
      "Site: 758\n",
      "Site: 829\n",
      "Site: 1282\n",
      "Site: 1650\n",
      "Site: 1770\n",
      "Site: 1800\n",
      "Site: 2077\n",
      "Site: 2899\n",
      "Site: 2974\n",
      "Site: 2985\n",
      "Site: 2993\n",
      "Site: 3071\n",
      "Site: 3083\n",
      "Site: 3188\n",
      "Site: 3529\n",
      "Site: 3547\n",
      "Site: 3593\n",
      "Site: 3851\n",
      "Site: 3868\n",
      "Site: 3904\n",
      "Site: 4639\n",
      "Total number of cuts: 29\n"
     ]
    }
   ],
   "source": [
    "from Bio import SeqIO\n",
    "from Bio.Restriction import RestrictionBatch\n",
    "\n",
    "# Load the plasmid sequence file\n",
    "plasmid_file = \"Desktop/pGL3.fa\"\n",
    "\n",
    "# Read the plasmid sequence using Biopython SeqIO\n",
    "plasmid_record = SeqIO.read(plasmid_file, \"fasta\")\n",
    "\n",
    "# Define the restriction enzymes and their recognition sites\n",
    "enzymes = {\n",
    "    \"ApoI\": \"RAATTY\",\n",
    "    \"BsaI\": \"GGTCTC\",\n",
    "    \"DpnI\": \"GATC\"\n",
    "}\n",
    "\n",
    "# Create a RestrictionBatch object with the enzymes\n",
    "batch = RestrictionBatch(enzymes)\n",
    "\n",
    "# Perform the restriction digest on the plasmid sequence\n",
    "digest_data = batch.search(plasmid_record.seq)\n",
    "\n",
    "# Display the digest results\n",
    "for enzyme, sites in digest_data.items():\n",
    "    print(f\"Results for {enzyme}:\")\n",
    "    for site in sites:\n",
    "        print(f\"Site: {site}\")\n",
    "\n",
    "# Total number of cuts\n",
    "total_cuts = sum(len(sites) for sites in digest_data.values())\n",
    "print(f\"Total number of cuts: {total_cuts}\")\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "d13b2c82-dc46-4a74-b8df-31b153cd1bb4",
   "metadata": {},
   "source": [
    "Restriction digest script - matches forward and reverse strands of all three enzyme sites.\n",
    "Code should include a regex sequence to match the three restriction sites on both strands of the DNA. This can be accomplished by 1) having both forward and reverse motifs in the regular expression, 2) searching forward strand, reverse complementing, and searching again. (The code uses re.finditer method here.)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "e2894464-882e-4db6-afb6-77cbf910aae8",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Number of matches for ApoI: 4\n",
      "Locations of matches for ApoI:\n",
      "Start: 573, End: 579\n",
      "Start: 1929, End: 1935\n",
      "Start: 4567, End: 4573\n",
      "Start: 4578, End: 4584\n",
      "Number of matches for BsaI: 0\n",
      "Locations of matches for BsaI:\n",
      "Number of matches for DpnI: 24\n",
      "Locations of matches for DpnI:\n",
      "Start: 119, End: 123\n",
      "Start: 126, End: 130\n",
      "Start: 751, End: 755\n",
      "Start: 768, End: 772\n",
      "Start: 839, End: 843\n",
      "Start: 1292, End: 1296\n",
      "Start: 1660, End: 1664\n",
      "Start: 1780, End: 1784\n",
      "Start: 1810, End: 1814\n",
      "Start: 2087, End: 2091\n",
      "Start: 2909, End: 2913\n",
      "Start: 2984, End: 2988\n",
      "Start: 2995, End: 2999\n",
      "Start: 3003, End: 3007\n",
      "Start: 3081, End: 3085\n",
      "Start: 3093, End: 3097\n",
      "Start: 3198, End: 3202\n",
      "Start: 3539, End: 3543\n",
      "Start: 3557, End: 3561\n",
      "Start: 3603, End: 3607\n",
      "Start: 3861, End: 3865\n",
      "Start: 3878, End: 3882\n",
      "Start: 3914, End: 3918\n",
      "Start: 4649, End: 4653\n"
     ]
    }
   ],
   "source": [
    "import re\n",
    "\n",
    "# Read the plasmid sequence from the file\n",
    "with open(\"pGL3.fa\", \"r\") as file:\n",
    "    plasmid_sequence = file.read().replace(\"\\n\", \"\")\n",
    "\n",
    "# Define the restriction enzymes and their recognition sites\n",
    "enzymes = {\n",
    "    \"ApoI\": \"RAATTY\",\n",
    "    \"BsaI\": \"GGTCTC\",\n",
    "    \"DpnI\": \"GATC\"\n",
    "}\n",
    "# Function to find all matches of a regex motif in the plasmid sequence\n",
    "def find_all_matches(regex_pattern):\n",
    "    matches = [(match.start(), match.end()) for match in re.finditer(regex_pattern, plasmid_sequence)]\n",
    "    return matches\n",
    "\n",
    "# Find and display all matches for each restriction enzyme using regex\n",
    "for enzyme, regex_pattern in enzyme_regex.items():\n",
    "    matches = find_all_matches(regex_pattern)\n",
    "    num_matches = len(matches)\n",
    "    print(f\"Number of matches for {enzyme}: {num_matches}\")\n",
    "    print(f\"Locations of matches for {enzyme}:\")\n",
    "    for match in matches:\n",
    "        print(f\"Start: {match[0]}, End: {match[1]}\")\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "682c04c0-cbde-49df-ba31-b5c83d6d1615",
   "metadata": {},
   "source": [
    "The code uses re.findall() method from the Python re module to find all matches of a regex motif in the plasmid sequence"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "a7823f8f-1cee-461d-8533-c47ce7b0dfcc",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Number of matches for ApoI: 4\n",
      "Locations of matches for ApoI:\n",
      "Start: 573, End: 579\n",
      "Start: 573, End: 579\n",
      "Start: 573, End: 579\n",
      "Start: 4578, End: 4584\n",
      "Number of matches for BsaI: 0\n",
      "Locations of matches for BsaI:\n",
      "Number of matches for DpnI: 24\n",
      "Locations of matches for DpnI:\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n",
      "Start: 119, End: 123\n"
     ]
    }
   ],
   "source": [
    "import re\n",
    "\n",
    "# Read the plasmid sequence from the file\n",
    "with open(\"pGL3.fa\", \"r\") as file:\n",
    "    plasmid_sequence = file.read().replace(\"\\n\", \"\")\n",
    "\n",
    "# Define the restriction enzymes and their recognition sites\n",
    "enzymes = {\n",
    "    \"ApoI\": \"RAATTY\",\n",
    "    \"BsaI\": \"GGTCTC\",\n",
    "    \"DpnI\": \"GATC\"\n",
    "}\n",
    "# Function to find all matches of a regex motif in the plasmid sequence\n",
    "def find_all_matches(regex_pattern):\n",
    "    return re.findall(regex_pattern, plasmid_sequence)\n",
    "\n",
    "# Find and display all matches for each restriction enzyme using findall\n",
    "for enzyme, regex_pattern in enzyme_regex.items():\n",
    "    matches = find_all_matches(regex_pattern)\n",
    "    num_matches = len(matches)\n",
    "    print(f\"Number of matches for {enzyme}: {num_matches}\")\n",
    "    print(f\"Locations of matches for {enzyme}:\")\n",
    "    for match in matches:\n",
    "        print(f\"Start: {plasmid_sequence.find(match)}, End: {plasmid_sequence.find(match) + len(match)}\")\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "977989ca-2a6d-4295-b971-fcc8352e98e2",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Start: 119, End: 123\n"
     ]
    }
   ],
   "source": [
    "print(f\"Start: {plasmid_sequence.find(match)}, End: {plasmid_sequence.find(match) + len(match)}\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "d1ac64e5-ab33-44aa-9d8a-0edde98dd9f9",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
