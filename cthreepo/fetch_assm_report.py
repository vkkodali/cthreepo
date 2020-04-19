import requests
import re
import xml.etree.ElementTree as ET

eutils = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
esearch = eutils + 'esearch.fcgi'
esummary = eutils + 'esummary.fcgi'

# accn = 'GCF_000001635.26'

def perform_esearch(accn):
    acc_pattern = re.compile('^GC[FA]_\\d{9}\\.\\d{1,2}$')
    accn = accn.upper().strip()
    assert acc_pattern.match(accn), 'Incorrect Assembly Accession'
    esearch_params = {
        'db' : 'assembly', 
        'term' : accn + '[Assembly Accession]',
        'usehistory' : 'y'}
    try:
        r = requests.get(esearch, params=esearch_params)
        res = ET.fromstring(r.text)
    except requests.exceptions.RequestException as e:
        raise SystemExit(e)
    return res 

def perform_esummary(res):
    assert res.find('./Count').text == '1', 'esearch count >1'
    esummary_params = {
        'db' : 'assembly',
        'query_key' : res.find('./QueryKey').text,
        'WebEnv' : res.find('./WebEnv').text }
    try:
        r = requests.get(esummary, params=esummary_params)
        res_out = ET.fromstring(r.text)
    except requests.exceptions.RequestException as e:
        raise SystemExit(e)
    return res_out 

def fetch_asm_rpt(res):
    assm_rpt_path = res.find('.//FtpPath_Assembly_rpt').text
    assm_rpt_path = assm_rpt_path.replace('ftp://', 'https://')
    try:
        r = requests.get(assm_rpt_path)
    except requests.exceptions.RequestException as e:
        raise SystemExit(e)
    return r.text.splitlines()


