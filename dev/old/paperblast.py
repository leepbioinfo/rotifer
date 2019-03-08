#! /usr/bin/env python3
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from bs4 import BeautifulSoup
import warnings
warnings.filterwarnings("ignore")
import os
import argparse
import sys
sys.path.insert(0, '/home/kaihami/mymodules')

import rotifer.core.cli as corecli

__version__ = '0.2'
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = corecli.parser(description ='PaperBlast')
    parser.add(rconfig=':cli.acc')
    args = parser.parse_args()
    return args




def send_form(acc):
    form = driver.find_element_by_xpath('/html/body/form/p[1]/textarea') # form area

    # print('## sending', acc)
    form.send_keys(acc) # input acc
    form.submit() # submit
    try:
        # print('## is form')
        ok = driver.find_element_by_xpath(form)
    except:
#        print('## ok')
        return 'ok'

def go_soup(html, idx):
    soup = BeautifulSoup(html, 'html') # parse using BS
    #all_p = soup.findAll('p')
    events = {'onmousedown', 'onmouseup', 'onclick'}
    p_tags = soup.find_all(lambda tag: any(attr in events for attr in tag.attrs.keys()))

    f = []
    for tag in p_tags:
        tag2 = tag
        tag = str(tag)
        if 'curated' in str(tag):
            locus = '.'
    #        pid = '.'
            organism = '.'
            link = '.'
            description = '.'
            # Header
            _ = tag.split("'")
            #pid = tag2['href'].split('/')[-1]
            locus = [x for x in _ if 'curated' in x][0].split(':')[-1]
            res = []
            res2 = []
            #res.append(pid)
            res.append(str(idx))
            res.append(locus)

            for q in tag2.find_next_siblings():
                if 'of query is similar to' in str(q):
                    _2 = q.text.split(',')
                    identity = str(int(_2[0].strip().split(' ')[0].replace('%',''))/100)
                    coverage = str(int(_2[1].strip().split(' ')[0].replace('%', ''))/100)
                    res.append(identity)
                    res.append(coverage)
                if '<i>' in str(q):
                    organism = str(q).strip().replace('<i>', '').replace('</i>', '')
                    res.append(organism)
        if 'pb:' in str(tag):
            link = tag2['href']
            link_text = tag2.text
            for t in tag2.find_next_siblings():
                if 'ul' in str(t) and 'small' not in str(t):
                    description = t.text
                    # locus, ID, COV, link, title, description
                    print('\t'.join(res + [link_text,description, link]))


def ipt_accs(accs):
    idx = 0
    for acc in accs:
        idx +=1
        driver.get('http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi') # paperblast site
        form_sender = True
        while form_sender:
            send = send_form(acc)
            if send == 'ok':
                form_sender = False
        html = driver.page_source # parse page
        go_soup(html, idx)

if __name__ == '__main__':
    try:
        driver = webdriver.PhantomJS(executable_path='/home/kaihami/projects/gen_scripts/paperblast/phantomjs-2.1.1-linux-x86_64/bin/phantomjs')
        args = parse_cli()
        accs = args.accession
        print('\t'.join('Index Locus Identity Coverage Title Description Link'.split()))
        ipt_accs(accs)
    except KeyboardInterrupt:
        sys.exit(0)
