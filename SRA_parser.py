##next page
##get links
import urllib.request
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.support.ui import Select
import time
import pandas as pd  # specify the url

gse = "GSE13750"
urlpage = "https://www.ncbi.nlm.nih.gov/bioproject/?term=" + gse
driver = webdriver.Firefox()
driver.get(urlpage)
links = []
dic = {}
dic[gse] = {}
dic[gse]["pubmed"] = driver.find_element_by_id("PubMed").get_attribute("href")
dic[gse]["desc"] = driver.find_element_by_css_selector(".Title > h3:nth-child(2)").text
results = driver.find_element_by_id("SRA Experiments").click()
driver.find_element_by_css_selector(
    "div.results_settings:nth-child(1) > ul:nth-child(1) > li:nth-child(2) > a:nth-child(1)").click()
driver.find_element_by_css_selector("input[type='radio'][value='200']").click()
# nr_pages = driver.find_element_by_id("pageno").get_attribute("last")
# print(nr_pages)
hasNextPage = True
next_page = driver.find_elements_by_css_selector("div.pagination:nth-child(1) > a:nth-child(4)")
links = links + [x.get_attribute("href") for x in driver.find_elements_by_xpath("//*[contains(@class, 'title')]/a")]
while (len(next_page) > 0):
    next_page[0].click()
    links = links + [x.get_attribute("href") for x in driver.find_elements_by_xpath("//*[contains(@class, 'title')]/a")]
    next_page = driver.find_elements_by_css_selector("div.pagination:nth-child(1) > a:nth-child(4)")

for link in links:
    dic[gse][link] = {}
    driver.get(link)
    sra = [x.text for x in
           driver.find_elements_by_xpath("//a[contains(@href, 'trace.ncbi.nlm.nih.gov/Traces/sra/?run=')]")]
    dic[gse][link]["sra"] = sra
    desc = driver.find_elements_by_css_selector(".details > b:nth-child(1)")[0].text
    dic[gse][link]["desc2"] = \
    driver.find_elements_by_css_selector("div.sra-full-data:nth-child(3) > span:nth-child(1)")[0].text.split("\n")[0]
    # dic[gse][link]["desc"] = desc.split(":")[2].split(";")[0]

# print(links)
