import forestry
from bokeh.io import show,output_file

nameList = ['a','b','c','d','e','f','g','h','i','j']
parentList = [None,'h','a','c','b','d','e','a','c',None]
weightList = [0,1,2,6,4,8,3,9,3,0]
descList = ['This is','an amazing','example','of how','I can','add','things to','a description','field in','bokeh!']
setA=[nameList,parentList,weightList,descList]

nameList = ['EB_23', 'EB_24', 'TE_86', 'YST_20', 'EB_33', 'CIS_32', 'TE_72', 'EC_75', 'TE_84', 'EB_22', 'EC_19', 'EB_87', 'EC_85', 'CIS_73', 'YST_88', 'EB_26', 'EB_25', 'TE_74']
parentList = ['TE_74', 'TE_74', 'EC_75', 'TE_72', 'EC_85', 'Precursor', 'TE_84', 'Precursor', 'TE_74', 'EC_85', 'YST_20', 'EC_85', 'EC_75', 'Precursor', 'EB_33', 'EC_85', 'EC_75', 'EC_85']
weightList = [401.0, 619.0, 367.0, 898.0, 429.0, 294.0, 582.0, 13.0, 432.0, 342.0, 654.0, 304.0, 158.0, 1694.0, 485.0, 320.0, 431.0, 274.0]
descList = [['+ LOH 1q'], [], [], [], [], [], ['+ SNV19', '+ SNV25', '+ SNV30'], ['+ SNV3', '+ SNV4', '+ SNV7', '+ SNV8', '+ SNV10', '+ SNV11', '+ SNV13', '+ SNV14', '+ SNV17', '+ SNV18', '+ SNV20', '+ SNV22', '+ SNV29', '+ LOH 4q', '+ LOH 14q', '+ LOH 15q', '+ LOH 22q'], ['+ SNV16', '+ SNV21'], [], ['+ SNV6'], [], [], [], [], [], [], []]
setB=[nameList,parentList,weightList,descList]

nameList = ['CIS29', 'NS48', 'BCIS28', 'NS30', 'FCIS27', 'NS75']
parentList = ['BCIS28', 'NS30', 'Precursor', 'CIS29', 'Precursor', 'NS48']
weightList = [1223.0, 846.0, 129.0, 653.0, 89.0, 709.0]
descList = [['+ SNV13', '+ SNV17'], [], ['+ SNV0', '+ SNV1', '+ SNV14', '+ SNV19', '+ SNV36', '+ SNV39', '+ SNV42', '+ LOH 5q'], ['+ SNV6', '+ SNV7', '+ SNV8', '+ SNV10', '+ SNV11', '+ SNV12', '+ SNV15', '+ SNV16', '+ SNV18', '+ SNV23', '+ SNV26', '+ SNV27', '+ SNV28', '+ SNV30', '+ SNV31', '+ SNV33', '+ SNV35', '+ SNV37', '+ SNV40', '+ SNV41', '+ LOH 1q', '+ LOH 9q', '+ LOH 22q'], ['+ SNV0', '+ SNV1', '+ SNV14', '+ SNV19', '+ SNV36', '+ SNV39', '+ SNV42', '+ LOH 4q', '+ LOH 5q'], ['+ SNV9', '+ SNV29']]
setC=[nameList,parentList,weightList,descList]

print [setA,setC,setB]
# optional:
output_file('trees.html', title='Make Forests Great Again!', mode='cdn', root_dir=None)
# required:
show(forestry.growForest([setA,setC,setB],['Simpele test','Hele domme boom','Minder domme boom']))