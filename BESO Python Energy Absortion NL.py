import math,customKernel
from abaqus import getInput,getInputs
from odbAccess import openOdb
from abaqusConstants import *
from timeit import default_timer as timer
## Função para execução da análise de elementos Finitos, cálculo do número de sensibilidade e da função objetivo
def FEA(Iter,Mdb,Xe,Ae,Nv,Fh,Elmts):
    try:
        Mdb.Job('Design_Job'+str(Iter),'Model-1').submit()
        Mdb.jobs['Design_Job'+str(Iter)].waitForCompletion()
    except AbaqusException, message:
        print "Error occured: ", message
        sys.exit(1)
    opdb = openOdb('Design_Job'+str(Iter)+'.odb')
    w = opdb.steps['Step-1'].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLWK'].data[-1][1]
    seng = opdb.steps['Step-1'].frames[-1].fieldOutputs['SENER'].values
    sensi_el = {}
    for en in seng: sensi_el[en.elementLabel] = en.data
    seng = opdb.steps['Step-1'].frames[-1].fieldOutputs['PENER'].values
    sensi_pl = {}
    for en in seng: sensi_pl[en.elementLabel] = en.data
    for en in seng:
        Ae[en.elementLabel] = (Xe[en.elementLabel]/len(Xe)) / Nv - (sensi_el[en.elementLabel] + sensi_pl[en.elementLabel]) / w
    obj = (w / (Nv * len(Elmts) * 0.25)) * 1000**2
    ft = {}
    fmax = []
    seng = opdb.steps['Step-1'].frames
    for en in seng:
        sing = en.fieldOutputs['RT'].values
        for yn in sing: ft[yn.nodeLabel] = yn.magnitude
        fmax.append(max(ft.values()))
    Fh.append(max(fmax))
    opdb.close()
    return obj
## Função para acessar as distâncias entre centros dos elementos calculadas separadamente pela função preFlt(Rmin,Elmts,Nds,Fm)
def preFm():
    mddbb = openMdb('Filtred_design.cae')
    fmm = mddbb.customData.History['fmd']
    mddbb.close()
    return fmm
## Função de cálculo das distâncias entre centros dos elementos
def preFlt(Rmin,Elmts,Nds,Fm):
    import numpy as np
    # Calculo das coordenadas do centro do elemento
    elm, c0 = np.zeros(len(Elmts)), np.zeros((len(Elmts),3))
    for i in range(len(elm)):
        elm[i] = Elmts[i].label
        nds = Elmts[i].connectivity
        for nd in nds: c0[i] = np.add(c0[i],np.divide(Nds[nd].coordinates,len(nds)))
    # Fator de peso
    for i in range(len(elm)):
        Fm[elm[i]] = [[],[]]
        for j in range(len(elm)):
            dis = np.sqrt(np.sum(np.power(np.subtract(c0[i],c0[j]),2)))
            if dis<Rmin: 
                Fm[elm[i]][0].append(elm[j])
                Fm[elm[i]][1].append(Rmin - dis)
        Fm[elm[i]][1] = np.divide(Fm[elm[i]][1],np.sum(Fm[elm[i]][1]))
## Função de aplicação do filtro
def fltAe(Ae,Fm):
    raw = Ae.copy()
    for el in Fm.keys():
        Ae[el] = 0.0
        for i in range(len(Fm[el][0])): Ae[el]+=raw[Fm[el][0][i]]*Fm[el][1][i]
## Função BESO
def BESO(Vf,Xe,Ae,Part,Elmts):
    lo, hi = min(Ae.values()), max(Ae.values())
    tv = Vf*len(Elmts)
    while math.fabs((hi-lo)/hi) > 1.0e-5:
        th = (lo+hi)/2.0
        for key in Xe.keys(): Xe[key] = 1.0 if Ae[key]<th else 0.0
        if sum(Xe.values())-tv>0: hi = th
        else: lo = th
    # Rotulando elementos como sólidos ou vazios
    vlb, slb = [], []
    for el in Elmts:
        if Xe[el.label] == 1.0: slb.append(el.label)
        else: vlb.append(el.label)
    for el in Elmts:
        if Xe[el.label] != 1.0:
            Ae[el.label] = 0.0
    # Atribuindo elementos sólidos e vazios a cada seção
    Part.SectionAssignment(Part.SetFromElementLabels('ss',slb),'sldSec')
    Part.SectionAssignment(Part.SetFromElementLabels('vs',vlb),'voidSec')
## ====== MAIN PROGRAM ======
if __name__ == '__main__':
    # Definir parâmetros e entradas
    pars = (('MaxForce:','1'), ('Rmin:', '1'), ('ER:', '0.02'))
    mf,rmin,ert = [float(k) if k!=None else 0 for k in getInputs(pars,dialogTitle='Parameters')]
    if mf<=0 or rmin<0 or ert<=0: sys.exit()
    fm = preFm()
    mddb = openMdb(getInput('Input CAE file:',default='Test.cae'))
    # Inicialização do Design
    part = mddb.models['Model-1'].parts['Part-1']
    elmts, nds = part.elements, part.nodes
    oh, vh, fh, tfea, ta, ch = [], [], [], [], [], []
    xe, ae, oae = {}, {}, {}
    for el in elmts: xe[el.label] = 1.0
    # Iterações de otimização
    change, iter, obj = 1, -1, 0
    inicio_1 = timer()
    while iter < 130:
        iter += 1
        if iter == 0:
            vh.append(sum(xe.values()) / len(xe))
            nv = vh[0]
            # Rodar função FEA
            inicio_2 = timer()
            oh.append(FEA(iter, mddb, xe, ae, nv, fh, elmts))
            final_2 = timer()
            difer_2 = final_2 - inicio_2
            tfea.append(difer_2)
            # Aplicação do filtro e cálculo da média de sensibilidade entre a iteração atual e a iteração anterior
            if rmin > 0: fltAe(ae, fm)
            if iter > 0: ae = dict([(k, (ae[k] + oae[k]) / 2.0) for k in ae.keys()])
            oae = ae.copy()
            # Otimização BESO
            if fh[-1] > mf:
                nv = vh[-1] * (1.0 - ert)
            else:
                nv = vh[-1] * (1.0 + ert)
            BESO(nv, xe, ae, part, elmts)
            a = mddb.models['Model-1'].rootAssembly
            region = a.instances['Part-1-1'].sets['vs']
            mddb.models['Model-1'].ModelChange(name='Int-2', createStepName='Step-1',
                                               region=region, regionType=ELEMENTS, activeInStep=False,
                                               includeStrain=False)
        else:
            vh.append(sum(xe.values()) / len(xe))
            nv = vh[-1]
            # Rodar função FEA
            inicio_2 = timer()
            oh.append(FEA(iter, mddb, xe, ae, nv, fh, elmts))
            final_2 = timer()
            difer_2 = final_2 - inicio_2
            tfea.append(difer_2)
            # Aplicação do filtro e cálculo da média de sensibilidade entre a iteração atual e a iteração anterior
            if rmin>0: fltAe(ae,fm)
            if iter > 0: ae=dict([(k,(ae[k]+oae[k])/2.0) for k in ae.keys()])
            oae = ae.copy()
            # Otimização BESO
            if fh[-1] > mf:
                nv = vh[-1] * (1.0 - ert)
            else:
                nv = vh[-1] * (1.0 + ert)
            BESO(nv, xe, ae, part, elmts)
        if iter>10: change=math.fabs((sum(oh[iter-4:iter+1])-sum(oh[iter-9:iter-4]))/sum(oh[iter-9:iter-4]))
        ch.append(change)
    final_1 = timer()
    difer_1 = final_1 - inicio_1
    ta.append(difer_1)
    # Salvando resultados
    mddb.customData.History = {'vol':vh,'obj':oh,'tempfea':tfea,'tempa':ta,'change':ch,'force':fh}
    mddb.saveAs('Final_design.cae')
