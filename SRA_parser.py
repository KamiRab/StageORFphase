##next page
##get links
import urllib.request
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.support.ui import Select
import time
import pandas as pd  # specify the url

from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException
from selenium.common.exceptions import NoSuchElementException
# urlpage = "http://sysbio.gzzoc.com/rpfdb/browse.html#@Scerevisiae@GSE100626"
# driver = webdriver.Firefox()
# try:
#     driver.get(urlpage)
# except TimeoutException:
#     time.sleep(5)
#     driver.get(urlpage)
#
# dic = {}
# tabs = driver.find_elements_by_xpath("//ul[@id = 'ScerevisiaeStudyNav']/li[contains(@class, 'tab')]/a")
# # name = driver.find_element_by_class_name("tab active").text
#
# # dic[name]["pmid"] = driver.find_element_by_css_selector("#Scerevisiae_GSE100626_Detail > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(4) > a:nth-child(1)").get_attribute("href"))
# for tab in tabs:
#     link = tab.get_attribute("href")
#     tab.click()
#     WebDriverWait(driver,10).until(EC.presence_of_element_located((By.ID,"ScerevisiaeDetail")))
#     name = tab.text
#     dic[name]={}
#     dic[name]["pmid"] = driver.find_element_by_xpath("/html/body/div[2]/div/div[2]/div[1]/div[22]/div/div[@class='tab-pane active']/table/tbody/tr[2]/td[4]/a").get_attribute("href")
#     dic[name]["accession"] = driver.find_element_by_xpath("/html/body/div[2]/div/div[2]/div[1]/div[22]/div/div[@class='tab-pane active']/table/tbody/tr[2]/td[2]/a").get_attribute("href")
#     # dic[name]["strain"] = driver.find_element_by_xpath("/html/body/div[2]/div/div[2]/div[1]/div[22]/div/div[@class='tab-pane active']/table/tbody/tr[3]/td[2]").text.split(":")[1]
#     strain1= driver.find_element_by_xpath("/html/body/div[2]/div/div[2]/div[1]/div[22]/div/div[@class='tab-pane active']/table/tbody/tr[3]/td[2]").text
#     strain2= strain1.split(":")
#     try:
#         strain3= strain2[1]
#     except IndexError:
#         strain3 = strain2[0]
#     dic[name]["strain"] = strain3
# for GSE in dic:
#     url = dic[GSE]["accession"]
#     try:
#         driver.get(url)
#     except TimeoutException:
#         time.sleep(5)
#         driver.get(url)
#     dic[GSE]["pubmed"] = driver.find_element_by_id("PubMed").get_attribute("href")
#     dic[GSE]["desc"] = driver.find_element_by_css_selector(".Title > h3:nth-child(2)").text
#     results = driver.find_element_by_id("SRA Experiments").click()
#     driver.find_element_by_css_selector(
#         "div.results_settings:nth-child(1) > ul:nth-child(1) > li:nth-child(2) > a:nth-child(1)").click()
#     WebDriverWait(driver,100000).until(EC.element_to_be_clickable((By.CSS_SELECTOR,"input[type='radio'][value='200']")))
#     driver.find_element_by_css_selector("input[type='radio'][value='200']").click()
#     # nr_pages = driver.find_element_by_id("pageno").get_attribute("last")
#     # print(nr_pages)
#     hasNextPage = True
#     next_page = driver.find_elements_by_css_selector("div.pagination:nth-child(1) > a:nth-child(4)")
#     links = []
#     links = links + [x.get_attribute("href") for x in driver.find_elements_by_xpath("//*[contains(@class, 'title')]/a")]
#     while (len(next_page) > 0):
#         next_page[0].click()
#         links = links + [x.get_attribute("href") for x in driver.find_elements_by_xpath("//*[contains(@class, 'title')]/a")]
#         next_page = driver.find_elements_by_css_selector("div.pagination:nth-child(1) > a:nth-child(4)")
#
#     for link in links:
#         dic[GSE][link] = {}
#         try:
#             driver.get(link)
#         except TimeoutException:
#             time.sleep(5)
#             driver.get(link)
#         sra = [x.text for x in
#                driver.find_elements_by_xpath("//a[contains(@href, 'trace.ncbi.nlm.nih.gov/Traces/sra/?run=')]")]
#         dic[GSE][link]["sra"] = sra
#         desc = driver.find_elements_by_css_selector(".details > b:nth-child(1)")[0].text
#         dic[GSE][link]["desc2"] = \
#         driver.find_elements_by_css_selector("div.sra-full-data:nth-child(3) > span:nth-child(1)")[0].text.split("\n")[0]
#         # dic[gse][link]["desc"] = desc.split(":")[2].split(";")[0]

# print(links)
with open("SRA_parser.tab","w") as parser:
    parser.write("Accession\tAccession link\tPMID\tStrain\tPubmed\tdesc\tdesc2\tSRA experiments link\tSRA\n")
    urlpage = "http://sysbio.gzzoc.com/rpfdb/browse.html#@Scerevisiae@GSE100626"
    driver = webdriver.Firefox()
    try:
        driver.get(urlpage)
    except TimeoutException:
        time.sleep(5)
        driver.get(urlpage)

    dic = {}
    tabs = driver.find_elements_by_xpath("//ul[@id = 'ScerevisiaeStudyNav']/li[contains(@class, 'tab')]/a")
    # name = driver.find_element_by_class_name("tab active").text

    # dic[name]["pmid"] = driver.find_element_by_css_selector("#Scerevisiae_GSE100626_Detail > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(4) > a:nth-child(1)").get_attribute("href"))
    for tab in tabs:
        link = tab.get_attribute("href")
        tab.click()
        # time.sleep(10)
        WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, "ScerevisiaeDetail")))
        # WebDriverWait(driver, 10).until(EC.presence_of_all_elements_located)
        name = tab.text

        dic[name] = {}
        dic[name]["pmid"] = driver.find_element_by_xpath(
            "/html/body/div[2]/div/div[2]/div[1]/div[22]/div/div[@id='Scerevisiae_{}_Detail']/table/tbody/tr[2]/td[4]/a".format(name)).get_attribute(
            "href")
        dic[name]["accession"] = driver.find_element_by_xpath(
            "/html/body/div[2]/div/div[2]/div[1]/div[22]/div/div[@id='Scerevisiae_{}_Detail']/table/tbody/tr[2]/td[2]/a".format(name)).get_attribute(
            "href")
        dic[name]["title"] = driver.find_element_by_xpath("/html/body/div[2]/div/div[2]/div[1]/div[22]/div/div[@id='Scerevisiae_{}_Detail']/table/tbody/tr[1]/td[2]/a".format(name)).text
        # dic[name]["strain"] = driver.find_element_by_xpath("/html/body/div[2]/div/div[2]/div[1]/div[22]/div/div[@class='tab-pane active']/table/tbody/tr[3]/td[2]").text.split(":")[1]
        strain1 = driver.find_element_by_xpath(
            "/html/body/div[2]/div/div[2]/div[1]/div[22]/div/div[@id='Scerevisiae_{}_Detail']/table/tbody/tr[3]/td[2]".format(name).format(name)).text
        strain2 = strain1.split(":")
        try:
            strain3 = strain2[1]
        except IndexError:
            strain3 = strain2[0]
        dic[name]["strain"] = strain3
    for GSE in dic:
        url = dic[GSE]["accession"]
        pmid = dic[GSE]["pmid"]
        strain = dic[GSE]["strain"]
        try:
            driver.get(url)
        except TimeoutException:
            time.sleep(5)
            driver.get(url)
        if driver.find_element_by_tag_name("body").get_attribute("id") == "ui-ncbiexternallink-5":
            try:
                url = driver.find_element_by_xpath(
                    "//div[@class='rprt']/div[2]/p/a[text()='{}']".format(dic[GSE]["title"])).get_attribute("href")
            except NoSuchElementException:
                # url = driver.find_element_by_xpath(
                #     "//div[@class='rprt']/div[2]/div[1]/p[text()='{}']".format(dic[GSE]["title"])).get_attribute("href")
                try:
                    url = driver.find_element_by_xpath("//div[@class='rprt']/div[2]/div[@class = 'supp']/p[@class = 'desc'][text()='{}']//ancestor::div[@class = 'rslt']//child::p/a".format(dic[GSE]["title"])).get_attribute("href")
                except:
                    url = "Not found"
                    continue
            try:
                driver.get(url)
            except TimeoutException:
                time.sleep(5)
                driver.get(url)
        elif driver.find_element_by_tag_name("body").get_attribute("id") == "ui-ncbiexternallink-4":
            continue
        try:
            pubmed = driver.find_element_by_id("PubMed").get_attribute("href")
        except NoSuchElementException:
            pubmed = "not found"
        try:
            desc = driver.find_element_by_css_selector(".Title > h3:nth-child(2)").text
        except NoSuchElementException:
            desc = driver.find_element_by_css_selector(".Title > h2:nth-child(1)").text
        results = driver.find_element_by_id("SRA Experiments").get_attribute("href")
        driver.get(results)
        try:
            driver.find_element_by_css_selector(
                "div.results_settings:nth-child(1) > ul:nth-child(1) > li:nth-child(2) > a:nth-child(1)").click()
            WebDriverWait(driver, 100000).until(
                EC.element_to_be_clickable((By.CSS_SELECTOR, "input[type='radio'][value='200']")))
            driver.find_element_by_css_selector("input[type='radio'][value='200']").click()
        except NoSuchElementException :
            print("One page ? " + results)
            # nr_pages = driver.find_element_by_id("pageno").get_attribute("last")
            # print(nr_pages)
        next_page = driver.find_elements_by_css_selector("div.pagination:nth-child(1) > a:nth-child(4)")
        links = []
        try:
            links = links + [x.get_attribute("href") for x in
                             driver.find_elements_by_xpath("//*[contains(@class, 'title')]/a")]
            while (len(next_page) > 0):
                next_page[0].click()
                links = links + [x.get_attribute("href") for x in
                                 driver.find_elements_by_xpath("//*[contains(@class, 'title')]/a")]
                next_page = driver.find_elements_by_css_selector("div.pagination:nth-child(1) > a:nth-child(4)")
        except NoSuchElementException:
            links = [url]


        for sra_link in links:
            try:
                driver.get(sra_link)
            except TimeoutException:
                time.sleep(5)
                driver.get(sra_link)
            sra = [x.text for x in
                   driver.find_elements_by_xpath("//a[contains(@href, 'trace.ncbi.nlm.nih.gov/Traces/sra/?run=')]")]
            try:
                desc2 = \
                driver.find_elements_by_css_selector("div.sra-full-data:nth-child(3) > span:nth-child(1)")[
                    0].text.split("\n")[0]
            except IndexError:
                desc2 = "to look"
            for i in sra:
                parser.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(GSE,url,pmid,strain,pubmed,desc,desc2,sra_link,i))

