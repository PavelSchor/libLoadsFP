import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
import sys

global ll_dynPressure,ll_cpName,ll_dpName,ll_normalsNmt, ll_normalsConsist; ll_workingOn='CellData'
ll_dynPressure=1.
ll_cpName='cp';ll_dpName='dp'
ll_normalsNmt=False;
ll_normalsConsist=False

class LLSettings(object):
    def __init__(self):
        pass

    def fromGlobal(self,names):
        for name in names:
            try:
                setattr(self, name, getattr(sys.modules[__name__],name))
                # print(getattr(sys.modules[__name__],name), __name__)
            except:
                print(f'Can not set {name} from global llSetting')
                pass

    def fromOther(self,other,names):
        for name in names:
            try:
                setattr(self, name, getattr(other,name))
            except:
                pass
    def toGlobal(self,names):
        for name in names:
            try:
                setattr(sys.modules[__name__], name, getattr(self,name))
            except:
                pass
    def toOther(self,other,names):
            try:
                setattr(other, name, getattr(self,name))
            except:
                pass

    def print(self):
        print(vars(self))

def llSet(name,value=None):
    if isinstance(name,str):
        setattr(sys.modules[__name__], name,value)
    if isinstance(name,dict):
        for i in name:
            setattr(sys.modules[__name__], i,name[i])

def llGetKw(info,argL=[]):
    kw=LLSettings();
    if info is not None:
        # print(vars(info))
        kw.fromOther(info, list(vars(info).keys()) )
    else:
        kw.fromGlobal(argL)
    return kw

def llLoadVtk(fname):
    r=vtk.vtkPolyDataReader()
    r.SetFileName(fname)
    r.Update()
    # r.GetOutputPort()
    return r

def llSaveVtk(inp,fname):
    w=vtk.vtkPolyDataWriter()
    w.SetFileName(fname)
    w.SetInputConnection(inp.GetOutputPort())
    w.Update()

def llCellNormals(inp):
    #print(ll_normalsConsist,ll_normalsNmt)
    f=vtk.vtkPolyDataNormals()
    f.SetInputConnection(inp.GetOutputPort())
    f.SetComputeCellNormals(True)
    f.Update()
    return f

def llCellSize(inp):
    f=vtk.vtkCellSizeFilter()
    f.SetInputConnection(inp.GetOutputPort())
    f.Update()
    return f


def llCellCenters(inp):
    f=vtk.vtkCellCenters()
    f.SetInputConnection(inp.GetOutputPort())
    f.Update()
    xi=dsa.WrapDataObject(inp.GetOutput())
    xf=dsa.WrapDataObject(f.GetOutput())
    xi.CellData.append(xf.Points,'CellCenters')
    inp.Update()
    return inp


def llP2F(inp,info=None):

    kw=llGetKw(info,['ll_dynPressure','ll_cpName','ll_dpName'])
    #print(locals())
    #print(sys.modules[__name__].ll_dynPressure)
    #print(hasattr(sys.modules[__name__],'ll_dynPressure' ))
    q=sys.modules[__name__].ll_dynPressure;  cpn=sys.modules[__name__].ll_cpName;  dpn=sys.modules[__name__].ll_dpName
    # q=1,cpn=None,dpn=None

    q=kw.ll_dynPressure
    cpn=kw.ll_cpName
    dpn=kw.ll_dpName
    kw=None

    inp1=llCellSize(inp)
    inp2=llCellNormals(inp1)

    x=dsa.WrapDataObject(inp.GetOutput())
    dk=x.CellData.keys()
    # if not cpn in dk:
    #     print(f'{cpn} not in input array')
    # if not dpn in dk:
    #     print(f'{dpn} not in input array')

    f=vtk.vtkArrayCalculator()
    f.SetAttributeTypeToCellData()
    f.SetInputConnection(inp2.GetOutputPort())
    f.AddScalarArrayName(cpn)
    f.AddScalarArrayName('Area')
    f.AddVectorArrayName('Normals')
    f.SetFunction(f"Normals*{cpn}*Area*{q}")
    f.SetResultArrayName(f"Force_{cpn}")
    f.Update()
    return f

def llClipByPlanes(inp,info=None,plOrigins=[],plNormals=[],insideOut=None):
    n=len(plOrigins)
    kw=llGetKw(info,['ll_dynPressure','ll_cpName','ll_dpName','ll_clipInsideOut'])
    q=kw.ll_dynPressure
    cpn=kw.ll_cpName

    if insideOut is None:
        flip=kw.ll_clipInsideOut
    else:
        flip=insideOut
    kw=None
    if n == 0:
        return inp

    else:

        cpd={}
        pld={}

        pld[0]=vtk.vtkPlane()
        pld[0].SetNormal(plNormals[0])
        pld[0].SetOrigin(plOrigins[0])

        cpd[0]=vtk.vtkClipDataSet()
        cpd[0].SetInputConnection(inp.GetOutputPort())
        cpd[0].SetClipFunction(pld[0])
        cpd[0].Update()
        if n==1:
            return cpd[0]
        for i in range(1,n):
            pld[i]=vtk.vtkPlane()
            pld[i].SetNormal(plNormals[i])
            pld[i].SetOrigin(plOrigins[i])

            cpd[i]=vtk.vtkClipDataSet()
            cpd[i].SetInputConnection(cpd[i-1].GetOutputPort())
            cpd[i].SetClipFunction(pld[i])
            cpd[i].SetInsideOut(flip[i])
            cpd[i].Update()
        gf=vtk.vtkGeometryFilter()
        gf.SetInputConnection(cpd[i].GetOutputPort())
        gf.Update()
        return gf

def llComputeForceArraysAtPoint(inp,info=None,pt=[0,0,0]):
    kw=llGetKw(info,['ll_dynPressure','ll_cpName','ll_dpName'])
    q=kw.ll_dynPressure
    cpn=kw.ll_cpName
    kw=None
    cf=llP2F(inp,info)
    cc=llCellCenters(cf)

    f=vtk.vtkArrayCalculator()
    f.SetAttributeTypeToCellData()
    f.SetInputConnection(cc.GetOutputPort())
    # f.AddCoordinateVectorVariable( "coords",  0, 1, 2 )
    f.AddVectorArrayName('CellCenters')
    f.SetFunction(f"CellCenters-({pt[0]}*iHat+{pt[1]}*jHat+{pt[2]}*kHat)")
    f.SetResultArrayName(f"Arm")
    f.Update()

    g=vtk.vtkArrayCalculator()
    g.SetAttributeTypeToCellData()
    g.SetInputConnection(f.GetOutputPort())
    # g.AddCoordinateVectorVariable( "coords",  0, 1, 2 )
    g.AddVectorArrayName(f"Force_{cpn}")
    g.AddVectorArrayName(f'Arm')
    g.SetFunction(f"cross(Arm,Force_{cpn})")
    g.SetResultArrayName(f"Moment_{cpn}")
    g.Update()
    return g

def llGetArea(inp,info=None):
    f=llCellSize(inp)
    x=dsa.WrapDataObject(f.GetOutput())
    return x.CellData[f'Area'].sum()


def llGetForceArraysAtPoint(inp,info=None,pt=[0,0,0]):
    kw=llGetKw(info,['ll_dynPressure','ll_cpName','ll_dpName'])
    q=kw.ll_dynPressure
    cpn=kw.ll_cpName
    kw=None
    g=llComputeForceArraysAtPoint(inp,info,pt)
    x=dsa.WrapDataObject(g.GetOutput())
    pts=x.Points
    force=x.CellData[f'Force_{cpn}']
    moment=x.CellData[f'Moment_{cpn}']
    return pts, force, moment

def llGetForceAtPoint(inp,info=None,pt=[0,0,0]):
    pts,force,moment=llGetForceArraysAtPoint(inp,info,pt)
    return pt, force.sum(axis=0), moment.sum(axis=0)

def llGetForcesAtPointsAsArrays(inp,info=None,pts=None):
    n=len(inp)
    if pts is None:
        pts=np.zeros( (n,3))
    force=np.zeros((n,3))
    moment=np.zeros((n,3))
    for i in range(0,n):
        force[i,:], moment[i,:]= llGetForceAtPoint(inp[i],info,pts[i])[1:3]
    return pts, force, moment

def llGetForcesAtPointsAsPD(inp,info=None,pts=None):
    n=len(inp)
    if pts is None:
        pts=np.zeros( (n,3))
    force=np.zeros((n,3))
    moment=np.zeros((n,3))
    for i in range(0,n):
        force[i,:], moment[i,:]= llGetForceAtPoint(inp[i],info,pts[i])[1:3]
    p=vtk.vtkPolyData()
    x=dsa.WrapDataObject(p)
    x.Points=pts
    x.PointData.append(force,'force')
    x.PointData.append(moment,'moment')
    c=vtk.vtkAppendPolyData()
    c.SetInputData(p)
    c.Update()
    return c

def llWriteObjsDictAsVTKs(inp,info=None,fname='clip'):
    for i in inp:
        try:
            llSaveVtk(llComputeForceArraysAtPoint(inp[i]) ,f'{fname}_{i}.vtk')
        except:
            print(f'llWriteObjsDictAsVTKs: can not write {fname}_{i}.vtk')

def llGetClipsBy2PlanesByPoints(inp,info=None, cutPts=np.array([[0,0,0],[1,0,0]]), normals=np.array([[1,0,0]]) ):
    n=len(cutPts)-1
    cuts={}
    if len(normals)==1:
        normals=np.tile(normals,(n+1,1))
    for i in range(0,n):
        cuts[i]=llClipByPlanes(inp,info,np.array( [ cutPts[i], cutPts[i+1] ] ), np.array( [normals[i],normals[i+1] ]) ,None )
    return cuts

def llGetForceByClipsBy2PlanesAsPD(inp,info=None, cutPts=np.array([[0,0,0],[1,0,0]]), normals=np.array([[1,0,0]]),refPts=None ):
    if refPts is None:
        refPts=cutPts[1::]
    cuts=llGetClipsBy2PlanesByPoints(inp,info,cutPts, normals )
    return llGetForcesAtPointsAsPD(cuts,info,refPts)

def llGetForceByClipsBy2PlanesAsArrays(inp,info=None, cutPts=np.array([[0,0,0],[1,0,0]]), normals=np.array([[1,0,0]]),refPts=None ):
    if refPts is None:
        refPts=cutPts[1::]
    cuts=llGetClipsBy2PlanesByPoints(inp,info,cutPts, normals )
    return llGetForcesAtPointsAsArrays(cuts,info,refPts)


def llWriteClipsBy2PlanesAsVtk(inp,info=None, cutPts=np.array([[0,0,0],[1,0,0]]), normals=np.array([[1,0,0]]),refPts=None ,fname='clip'):
    cuts=llGetClipsBy2PlanesByPoints(inp,info,cutPts, normals )
    llWriteObjsDictAsVTKs(cuts,info,fname)

def llArrays2VTK(fname,coords,arrays):
    if not isinstance(arrays,dict):
        print(f'llArrays2VTK: arrays must be dict, not {type(arrays)}')

    p=vtk.vtkPolyData()
    x=dsa.WrapDataObject(p)
    x.Points=coords
    for i in arrays:
        x.PointData.append(arrays[i],i)
    c=vtk.vtkAppendPolyData()
    c.SetInputData(p)
    c.Update()
    llSaveVtk(c,fname)

def llGetAeroCoeffsRefArea(inp1,inp2,info=None, cutPts=np.array([[0,0,0],[1,0,0]]), normals=np.array([[1,0,0]]),refPts=None ):
    kw=llGetKw(info,['ll_dynPressure'])
    q=kw.ll_dynPressure
    kw=None
    if refPts is None:
        refPts=cutPts[1::]
    forces=llGetForceByClipsBy2PlanesAsArrays(inp1,info,cutPts,normals,refPts)
    cuts=llGetClipsBy2PlanesByPoints(inp2,info,cutPts, normals )
    areas=np.zeros((len(cuts),1))
    for i in range(0,len(cuts)):
        areas[i]=llGetArea(cuts[i])
    dist=np.linalg.norm(np.diff(cutPts,axis=0),axis=1)
    chords=areas/dist[None].T
    return forces[0], forces[1]/areas/q, forces[2]/areas/chords/q


