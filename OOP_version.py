import requests
import statistics
import pandas as pd
import numpy as np
import datetime 
from datetime import date
from dateutil.relativedelta import relativedelta

stock_list=["LKQ", "GPC", "PRTS", "AZO"]


    
url_get_financials = "https://apidojo-yahoo-finance-v1.p.rapidapi.com/stock/v2/get-financials"
url_get_details = "https://apidojo-yahoo-finance-v1.p.rapidapi.com/stock/get-detail"

headers = {'x-rapidapi-host': "apidojo-yahoo-finance-v1.p.rapidapi.com",
           'x-rapidapi-key': "181b5a206amsh0eee70c6e4ec573p113ea2jsnc87204bb7710"} 
my_bs_items=["intangibleAssets",
"capitalSurplus",
"totalLiab",
"totalStockholderEquity",
"otherCurrentLiab",
"totalAssets",
"endDate",
"commonStock",
"otherCurrentAssets",
"retainedEarnings",
"otherLiab",
"treasuryStock",
"otherAssets",
"cash",
"longTermDebt",
"shortTermInvestments",
"totalCurrentLiabilities",
"shortLongTermDebt",
"propertyPlantEquipment",
"totalCurrentAssets",
"netTangibleAssets",
"netReceivables"]

my_cf_items=["changeToLiabilities",
"totalCashflowsFromInvestingActivities",
"netBorrowings",
"totalCashFromFinancingActivities",
"changeToOperatingActivities",
"issuanceOfStock",
"netIncome",
"changeInCash",
"endDate",
"repurchaseOfStock",
"totalCashFromOperatingActivities",
"depreciation",
"changeToInventory",
"changeToAccountReceivables",
"otherCashflowsFromFinancingActivities",
"maxAge",
"changeToNetincome",
"capitalExpenditures"]

my_is_items=["researchDevelopment",
"effectOfAccountingCharges",
"incomeBeforeTax",
"minorityInterest",
"netIncome",
"sellingGeneralAdministrative",
"grossProfit",
"ebit",
"endDate",
"operatingIncome",
"otherOperatingExpenses",
"interestExpense",
"extraordinaryItems",
"nonRecurring",
"otherItems",
"incomeTaxExpense",
"totalRevenue",
"totalOperatingExpenses",
"costOfRevenue",
"totalOtherIncomeExpenseNet",
"maxAge",
"discontinuedOperations",
"netIncomeFromContinuingOps",
"netIncomeApplicableToCommonShares"]

sector_list={'Basic Materials':'XLB',
             'Communication Services':'XLC',
             'Consumer Cyclical':'XLY',
             'Consumer Defensive':'XLP',
             'Energy':'XLE',
             'Financial Services':'XLF',
             'Healthcare':'XLV', 
             'Industrials':'XLI',
             'Real Estate':'XLRE',
             'Technology':'XLK',
             'Utilities':'XLU'}

risk_free_rate=0.066
market_return=0.10
g = 0.025 

class Stock:


    def __init__(self,stock):
        self.stock = stock
    
    def run_api(self,url):
        querystring = {"region":"US","symbol": self.stock}
        response=requests.request("GET", url, headers=headers, params=querystring)
        response=response.json()
        return response
    
    def get_sector(self):
        sector=str(self['summaryProfile']['sector'].encode("utf-8").decode('unicode-escape').strip())
        return sector
    
    def shares_outstanding(self,stock):
        return self[stock]["defaultKeyStatistics"]["sharesOutstanding"]["raw"]
    

class FinancialStatements(Stock):

    
    def get_bs(self):
        return self['balanceSheetHistory']['balanceSheetStatements']
    
    def get_is(self):
        return self['incomeStatementHistory']['incomeStatementHistory']
         
    def get_cf(self):
        return self['cashflowStatementHistory']['cashflowStatements']  

class DcfAnalysis(Stock):
    
    def if_historical_odd(self,growth_type):
        for key in self:
            if self[key]<=0:
                print("For "+key+" :")
                print("The observed "+growth_type+ " is not positive")
                print("Historical growth rate I observed: ")
                print(self[key])
                self[key]=float(input("Please enter a 5 Year CAGR for projections:"))
            elif self[key]>=0.1 :
                print("For "+key+" :")
                print("The observed "+growth_type+ " is above 10%")
                print("Historical growth rate I observed: ")
                print(self[key])
                self[key]=float(input("Please enter a 5 Year CAGR for projections:"))
            else:
                self[key]=self[key]
        return self

    def fetch_row(self,data):
        fetched_row={}
        if isinstance(data, list)==True:
            for item in data:
                fetched_row[item]={k:[i[item]['raw'] for i in v] for (k,v) in self.items()}
        else:
            fetched_row[data]={k:[i[data]['raw'] for i in v] for (k,v) in self.items()}
        return fetched_row
     
    def yoy_growth_rates(self):
        return {key:{k:[(v[i] - v[i+1]) / float(v[i+1]) for i in range(0,len(v)-1)] 
        for (k,v) in value.items()} for (key,value) in self.items()}

    
    def mean_gr(self):
        return {key:{k:np.mean(v) for (k,v) in value.items()} for (key,value) in self.items()}

    def dict_in_dict_percentage(stockList,input_1,input_2):
       assumptions={}
       for stock in stockList:
           assumptions[stock]=[input_1[stock][i]/input_2[stock][i] for i in range(0,4)]
       return assumptions
    
    def del_extremes(self):
        output = [x for x in self if abs(x - np.mean(self)) < np.std(self) * 2]
        return output
        
    def second_highest(self):
        output = [x for x in self if x < max(self)]
        second_max= max(output)
        return second_max
    
    def second_lowest(self):
       output = [x for x in self if min(self)<x]
       second_min= min(output)
       return second_min
   
    def perpetuity_growth_dcf(self,periods,terminal_period):
        terminal_value = {k:v[-1] * (1 + g) / (wacc[k] - g) for (k,v) in self.items()}
        discount_rates_tv ={k:(1 / (1 + wacc[k])) ** v for (k,v) in terminal_period.items()}
        discount_rates = {k:[(1 / (1 + wacc[k])) ** i for i in v] for (k,v) in periods.items()}
        discounted_free_cash_flow={k:[v[i]* discount_rates[k][i] for i in range(1,5)] for (k,v) in self.items()}
        dcf_value = {k:sum(v) for (k,v) in discounted_free_cash_flow.items()}
        dcf_value = {k:v+terminal_value[k] * discount_rates_tv[k] for (k,v) in dcf_value.items()}
        return dcf_value
    
    def ebitda_exit_dcf(self,periods,ebitda,terminal_period):
        discount_rates_ebitda = {k:[(1 / (1 + wacc[k])) ** i for i in v] for (k,v) in periods.items()}
        discount_rates_tv ={k:(1 / (1 + wacc[k])) ** v for (k,v) in terminal_period.items()}
        terminal_ebitda = {k:v[4] for (k,v) in ebitda.items()}
        discounted_free_cash_flow_ebitda={k:[v[i]* discount_rates_ebitda[k][i] for i in range(1,5)] for (k,v) in self.items()}
        dcf_value_ebitda = {k:sum(v) for (k,v) in discounted_free_cash_flow_ebitda.items()}
        dcf_value_ebitda = {k:dcf_value_ebitda[k]+v*11*discount_rates_tv[k] for (k,v) in terminal_ebitda.items()}
        return dcf_value_ebitda

    def equity_value_dcf(self):
        equity_value= {k:v-total_debt[k]-bs_fetched_rows["cash"][k][0]-bs_fetched_rows["shortTermInvestments"][k][0] for (k,v) in self.items()}
        implied_share_price_dcf = {k:v/shares_outstanding[k] for (k,v) in equity_value.items()}
        print("The implied share price based on DCF: ")
        print(implied_share_price_dcf)
        return implied_share_price_dcf

    
def check_items_exist(data,item):
    for i in range(0,len(data)):
        data[i][item]=data[i].setdefault(item,0)
        if data[i][item]==0 or data[i][item]=={}:
            data[i][item]={'raw':0}
        else:
            item=item
    return(data)


get_financials={k:Stock(k).run_api(url_get_financials) for k in stock_list}

get_details={k:Stock(k).run_api(url_get_details) for k in stock_list}    
    
sector = {k:Stock.get_sector(v) for (k,v) in get_details.items()}

balance_sheet={k:FinancialStatements.get_bs(v) for (k,v) in get_financials.items()}

income_statement={k:FinancialStatements.get_is(v) for (k,v) in get_financials.items()}    

cash_flow={k:FinancialStatements.get_cf(v) for (k,v) in get_financials.items()}  
    
for i in my_bs_items:
    balance_sheet={k:check_items_exist(v,i) for (k,v) in balance_sheet.items()}  
    
for i in my_cf_items:
    cash_flow={k:check_items_exist(v,i) for (k,v) in cash_flow.items()}  

for i in my_is_items:
    income_statement={k:check_items_exist(v,i) for (k,v) in income_statement.items()}  
        


income_statement_items=['incomeTaxExpense', 'incomeBeforeTax','ebit']
cash_flows_items=['depreciation','capitalExpenditures']
balance_sheet_items=['totalCurrentAssets','cash','shortTermInvestments','longTermDebt', 'totalCurrentLiabilities', 'shortLongTermDebt']
rev_fetched_row=DcfAnalysis.fetch_row(income_statement,'totalRevenue')

is_fetched_rows= DcfAnalysis.fetch_row(income_statement,income_statement_items)
cf_fetched_rows= DcfAnalysis.fetch_row(cash_flow,cash_flows_items)
bs_fetched_rows= DcfAnalysis.fetch_row(balance_sheet,balance_sheet_items)

hist_yoy_gr=DcfAnalysis.yoy_growth_rates(rev_fetched_row)
hist_rev_gr_prof=DcfAnalysis.mean_gr(hist_yoy_gr)
{k:DcfAnalysis.if_historical_odd(v, "revenue growth rate") for (k,v) in hist_rev_gr_prof.items()} 

pro_forma_taxes=DcfAnalysis.dict_in_dict_percentage(stock_list,
                                    is_fetched_rows['incomeTaxExpense'],
                                    is_fetched_rows['incomeBeforeTax'])

pro_forma_taxes={k:np.mean(v) for (k,v) in pro_forma_taxes.items()}
pro_forma_taxes={k:v if v>0.15 and v<0.4 else 0.26 for (k,v) in pro_forma_taxes.items()}

dep_amor_percentage=DcfAnalysis.dict_in_dict_percentage(stock_list,
                                    cf_fetched_rows['depreciation'],
                                    rev_fetched_row['totalRevenue'])

dep_amor_percentage={k:DcfAnalysis.second_highest(v) for (k,v) in dep_amor_percentage.items()}

ebit=DcfAnalysis.dict_in_dict_percentage(stock_list,
                                    is_fetched_rows['ebit'],
                                    rev_fetched_row['totalRevenue'])

ebit_percentage={k:np.mean(DcfAnalysis.del_extremes(v)) for (k,v) in ebit.items()}


capex_percentage=DcfAnalysis.dict_in_dict_percentage(stock_list,
                                    cf_fetched_rows['capitalExpenditures'],
                                    rev_fetched_row['totalRevenue'])

capex_percentage={k:DcfAnalysis.second_highest(v) for (k,v) in capex_percentage.items()}

ebitda={stock:[is_fetched_rows['ebit'][stock][i]+cf_fetched_rows['depreciation'][stock][i] for i in range(0,4)] 
for stock in stock_list}

ebitda_percentage=DcfAnalysis.dict_in_dict_percentage(stock_list,
                                   ebitda,
                                    rev_fetched_row['totalRevenue'])

ebitda_percentage={k: DcfAnalysis.second_highest(v) for (k,v) in ebitda_percentage.items()}

net_working_cap={stock:[bs_fetched_rows['totalCurrentAssets'][stock][i]-
                        bs_fetched_rows['cash'][stock][i]-
                        bs_fetched_rows['shortTermInvestments'][stock][i]-
                        bs_fetched_rows['totalCurrentLiabilities'][stock][i]+
                        bs_fetched_rows['shortLongTermDebt'][stock][i]
                        for i in range(0,4)]  for stock in stock_list}

delta_net_working_cap={stock:[net_working_cap[stock][i+1]-
                             net_working_cap[stock][i] for i in range(0,3)]for stock in stock_list}

nwc_percentage=DcfAnalysis.dict_in_dict_percentage(stock_list,
                                   net_working_cap,
                                    rev_fetched_row['totalRevenue'])

nwc_percentage={k:v[0] for (k,v) in nwc_percentage.items()}

dates_historicals={k:[i['endDate']['fmt'] for i in v] for (k,v) in balance_sheet.items()}
mid_annual_filing_date= {k:datetime.datetime.strptime(v[0], '%Y-%m-%d').date()+relativedelta(months=-6)
for (k,v) in dates_historicals.items()}

filing_projections={k:[v+relativedelta(years=i) for i in range(0,6)] 
for (k,v) in mid_annual_filing_date.items()}

discount_periods={k:[abs(i-date.today()).days/float(365) for i in v[1:]] for (k,v) in filing_projections.items()}
sales = {k:pd.Series(index=v) for (k,v) in filing_projections.items()}

filing_projections_terminal={k:datetime.datetime.strptime(v[0], '%Y-%m-%d').date()+relativedelta(years=5) for (k,v) in dates_historicals.items()}
terminal_discount_periods={k:abs(v-date.today()).days/float(365)for (k,v) in filing_projections_terminal.items()}

for stock in stock_list:
    for i in range(1,6):
        sales[stock][0]=rev_fetched_row['totalRevenue'][stock][0]
        sales[stock][i]=sales[stock][i-1]*(1+hist_rev_gr_prof['totalRevenue'][stock]) 


capex_proj= {k:v* capex_percentage[k] for (k,v) in sales.items()}
ebitda_proj= {k:v*ebitda_percentage[k] for (k,v) in sales.items()}
ebit_proj= {k:v*ebit_percentage[k] for (k,v) in sales.items()}
nwc_proj= {k:v*nwc_percentage[k] for (k,v) in sales.items()}

for k in nwc_proj.keys():
    nwc_proj[k][0]=net_working_cap[k][0]
    
dep_amor_proj= {k:v*dep_amor_percentage[k] for (k,v) in sales.items()}
tax_proj= {k:v*pro_forma_taxes[k] for (k,v) in ebit_proj.items()}
   
   
delta_nwc_proj={k:[v[i-1]-v[i] for i in range(1,6)] for (k,v) in nwc_proj.items()}
     
fcf_proj={stock:[ebit_proj[stock][i]-
     tax_proj[stock][i]+
     capex_proj[stock][i]+
     dep_amor_proj[stock][i] for i in range(0,6)]for stock in stock_list}

fcf_proj={k:v[1:] for (k,v) in fcf_proj.items()}
fcf_proj={k:[v[i]+delta_nwc_proj[k][i] for i in range (0,5)] for (k,v) in fcf_proj.items()}

    
risk_free_rate= 0.066
market_return=0.10

shares_outstanding={k:v["defaultKeyStatistics"]["sharesOutstanding"]["raw"] for 
                    (k,v) in get_details.items()}

stock_beta={k:v["summaryDetail"]["beta"]["raw"] if v["summaryDetail"]["beta"]!={} else 0 
            for (k,v) in get_details.items()}

equity_value={k:v['price']['regularMarketPrice']["raw"]*shares_outstanding[k] for 
              (k,v) in get_financials.items()}

total_debt={k:v['financialData']['totalDebt']["raw"] if v["financialData"]["totalDebt"]!={} else 0 for (k,v) in get_details.items()}

return_on_equity={k:(0 if v==0 else risk_free_rate+v*(market_return-risk_free_rate)) for 
                  (k,v) in stock_beta.items()}

cost_of_debt={k:abs(v[0]["interestExpense"]["raw"])/total_debt[k] if 
              total_debt[k]!=0 else 0 for (k,v) in income_statement.items()}

wacc={k:0.10 if return_on_equity[k]==0 else
      risk_free_rate+v*(market_return-risk_free_rate) if v<0 else
      (equity_value[k]/(equity_value[k]+v))*return_on_equity[k]+(v/(equity_value[k]+v))*cost_of_debt[k]*(1-0.35) 
      for (k,v) in total_debt.items()}

g = 0.025
                    
# Perpetuity Growth DCF valuation
dcf_perpetuity=[]
dcf_ebitda_exit=[]

for k in sector.keys():
    if sector[k] == "Communication Services":
        dcf_perpetuity.append(k)
    elif sector[k] == "Technology":
        dcf_perpetuity.append(k)
    else:
        dcf_ebitda_exit.append(k)


dcf_perpetuity_growth_values=DcfAnalysis.perpetuity_growth_dcf(fcf_proj,discount_periods,terminal_discount_periods)
dcf_ebitda_exit_values=DcfAnalysis.ebitda_exit_dcf(fcf_proj,discount_periods,ebitda_proj,terminal_discount_periods)

perpetuity_growth_equity_value = DcfAnalysis.equity_value_dcf(dcf_perpetuity_growth_values)
ebitda_exit_equity_value = DcfAnalysis.equity_value_dcf(dcf_ebitda_exit_values)












