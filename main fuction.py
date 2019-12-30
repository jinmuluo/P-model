import Pmodel


a = Pmodel.pypmodel(tc=20, vpd=1000, co2=400, fapar=1, ppfd=300, elv=0, method_jmaxlim="smith19", c4='True',
                    do_soilmstress='True')
print(a)