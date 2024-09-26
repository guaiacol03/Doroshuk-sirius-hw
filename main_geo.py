import requests
from bs4 import BeautifulSoup
import bs4
from urllib.parse import urlparse, parse_qs
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str,
                    help='Card URL')
args = parser.parse_args()

url = args.input
try:
    res = requests.get(url)
except requests.exceptions.ConnectionError:
    raise Exception('URL invalid or unreachable')

if res.status_code != 200:
    raise Exception('Request rejected with ' + str(res.status_code))

soup = BeautifulSoup(res.text, 'html.parser')

# the only element with a consistent id is a search query
accession = soup.find("input", {"id": "geo_acc"})["value"] # get accession number
# table always contains header with accession number as id
strong_header = soup.find("strong", {"id": accession})
# table itself will have cellpadding attribute set
table = strong_header
while True:
    table = table.parent
    if table.has_attr("cellpadding"):
        break
# children of this table will be data columns
table_nodes = table.find_all(recursive=False)
table_dict = {}
for i, child in enumerate(table_nodes):
    k = child.find_all("td") # text fields are TR elements
    if len(k) != 2: # data fields contain exactly two
        continue
    table_dict[k[0].text] = (i, k[1].text) # first field is name and second - value

print("For experiment " + accession)
print("\tSpecies: " + table_dict["Organism"][1])
print("\tExperiment processing: " + table_dict["Treatment protocol"][1])
print("\tExperiment type: " + table_dict["Library strategy"][1])
print("\tGenotype: " + table_dict["Extracted molecule"][1])

sra_e = table_nodes[table_dict["SRA"][0]]
sra = sra_e.find("a")["href"]
print("\tSRA link: " + sra)

raw_e = table.find(text=re.compile(r"Series.*")).parent.parent.find_all("a")
print("\tRaw data links:")
for raw in raw_e:
    rval = raw["href"]
    print("\t\thttps://www.ncbi.nlm.nih.gov" + rval)