# import pandas as pd

# xls = pd.ExcelFile('enzymes.xlsx')

# df2 = pd.read_excel(xls, 'Sheet2')

# final_dict = {}

# for x in range(76):
#     y = {
#         df2['enzyme'][x]: {
#             "1bp": df2['1bp'][x],
#             "2bp": df2['2bp'][x],
#             "3bp": df2['3bp'][x],
#             "4bp": df2['4bp'][x],
#             "5bp": df2['5bp'][x],
#         }
#     }
#     final_dict.update(y)

# print(final_dict)


x = {
   "AccI":{
      "1bp":"-",
      "2bp":"-",
      "3bp":"-",
      "4bp":"-",
      "5bp":"-"
   },
   "AciI":{
      "1bp":"-",
      "2bp":"+",
      "3bp":"+",
      "4bp":"++",
      "5bp":"+++"
   },
   "AgeI-HF":{
      "1bp":"++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "AleI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "AluI":{
      "1bp":"-",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "ApaI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "AscI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "AvrII":{
      "1bp":"++",
      "2bp":"++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BamHI":{
      "1bp":"+",
      "2bp":"++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BamHI-HF":{
      "1bp":"+",
      "2bp":"+",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BbsI-HF":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BclI-HF":{
      "1bp":"-",
      "2bp":"-",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BglII":{
      "1bp":"++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BmtI-HF":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BsaI-HF":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BsiWI":{
      "1bp":"++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BsiWI-HF":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BsmBI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BssHII":{
      "1bp":"+",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "BstZ17I-HF":{
      "1bp":"+",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "ClaI":{
      "1bp":"-",
      "2bp":"-",
      "3bp":"+",
      "4bp":"+++",
      "5bp":"+++"
   },
   "DdeI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "DpnI":{
      "1bp":"-",
      "2bp":"++",
      "3bp":"++",
      "4bp":"nt",
      "5bp":"nt"
   },
   "DraIII-HF":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "EagI-HF":{
      "1bp":"+",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "EcoRI":{
      "1bp":"+",
      "2bp":"+",
      "3bp":"++",
      "4bp":"++",
      "5bp":"+++"
   },
   "EcoRI-HF":{
      "1bp":"+",
      "2bp":"+",
      "3bp":"++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "EcoRV":{
      "1bp":"++",
      "2bp":"++",
      "3bp":"++",
      "4bp":"++",
      "5bp":"+++"
   },
   "EcoRV-HF":{
      "1bp":"+",
      "2bp":"++",
      "3bp":"++",
      "4bp":"++",
      "5bp":"+++"
   },
   "Esp3I":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "FseI":{
      "1bp":"+",
      "2bp":"++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "HindIII":{
      "1bp":"-",
      "2bp":"+",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "HindIII-HF":{
      "1bp":"-",
      "2bp":"+",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "HpaI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "KpnI-HF":{
      "1bp":"+",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "MfeI-HF":{
      "1bp":"+",
      "2bp":"++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "MseI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "NcoI":{
      "1bp":"-",
      "2bp":"++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "NcoI-HF":{
      "1bp":"+",
      "2bp":"++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "NdeI":{
      "1bp":"+",
      "2bp":"+",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "NheI-HF":{
      "1bp":"++",
      "2bp":"++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "NlaIII":{
      "1bp":"++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "NotI":{
      "1bp":"++",
      "2bp":"++",
      "3bp":"++",
      "4bp":"++",
      "5bp":"++"
   },
   "NotI-HF":{
      "1bp":"++",
      "2bp":"++",
      "3bp":"++",
      "4bp":"++",
      "5bp":"++"
   },
   "NsiI":{
      "1bp":"+",
      "2bp":"+",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "NspI":{
      "1bp":"-",
      "2bp":"-",
      "3bp":"+",
      "4bp":"+",
      "5bp":"+++"
   },
   "PmeI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "PsiI":{
      "1bp":"+",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "PstI":{
      "1bp":"+",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "PstI-HF":{
      "1bp":"++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "PvuI-HF":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "PvuII":{
      "1bp":"++",
      "2bp":"++",
      "3bp":"++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "PvuII-HF":{
      "1bp":"-",
      "2bp":"++",
      "3bp":"++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "RsaI":{
      "1bp":"+",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SacI-HF":{
      "1bp":"-",
      "2bp":"+",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SacII":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SalI":{
      "1bp":"-",
      "2bp":"++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SalI-HF":{
      "1bp":"-",
      "2bp":"++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SapI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "Sau3AI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SbfI-HF":{
      "1bp":"++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "ScaI-HF":{
      "1bp":"+",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SfiI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SmaI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SpeI-HF":{
      "1bp":"+",
      "2bp":"++",
      "3bp":"++",
      "4bp":"++",
      "5bp":"++"
   },
   "SphI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SphI-HF":{
      "1bp":"++",
      "2bp":"++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "SspI-HF":{
      "1bp":"+",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "StuI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "StyI-HF":{
      "1bp":"+",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "XbaI":{
      "1bp":"++",
      "2bp":"++",
      "3bp":"++",
      "4bp":"++",
      "5bp":"++"
   },
   "XhoI":{
      "1bp":"++",
      "2bp":"++",
      "3bp":"++",
      "4bp":"+++",
      "5bp":"+++"
   },
   "XmaI":{
      "1bp":"+++",
      "2bp":"+++",
      "3bp":"+++",
      "4bp":"+++",
      "5bp":"+++"
   }
}