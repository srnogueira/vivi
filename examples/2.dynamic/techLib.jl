se_in = [
    Resource("steam",318.35,"kg/h",[0]),
    Resource("power",888.61+1.26,"kW",[0])
    ]

se_out = [
    Resource("h2",26.00,"kg/h",[0]),   # dry basis
    Resource("water",71.799,"kg/h",[0]),
    Resource("o2",206.38,"kg/h",[0])
]

se_heat = [
    HeatStruct(120.61,375.9,1023.2),
    HeatStruct(95.46,1023.2,473.2),
    HeatStruct(9.03,480.6,313.2),
    HeatStruct(61.98,480.6,1023.2),
    HeatStruct(105.61,1023.2,343.5),
    HeatStruct(55.68,343.5,313.2),
]

soec = Tech("soec",se_in,se_out,se_heat)

# Steam boiler
boiler_in = [
    Resource("power",0.017,"kW",[0]),
    Resource("water",3600,"kg/h",[0])
    ]
boiler_out = [Resource("steam",3600,"kg/h",[0])]
                
boiler_heat = [
    HeatStruct(326.1,298,376),
    HeatStruct(2249,376,377)
    ]

boiler = Tech("boiler",boiler_in,boiler_out,boiler_heat)

# Eletric heating
eh_in = [Resource("power",1000,"kW",[0])]
eh_heat = [HeatStruct(1000,1273,1272)]
eh = Tech("EH",eh_in,[],eh_heat)

# Cooling tower
ct_in = [Resource("power",14.3,"kW",[0])]
ct_heat = [HeatStruct(1000,298,318)]
ct = Tech("CT",ct_in,[],ct_heat)