using Luxor, Colors, HDF5
# Plotting mesh from JuliaFEM h5-file
# Currently deformed mesh plots only 1 mode of a modal solution.

function triangle(points, mesh_col, br)
	#Plotting a coloured triangle from points.
    jet = RGB{Float64}[
          RGB(clamp(min(4x - 1.5, -4x + 4.5) ,0.0,1.0*br),
              clamp(min(4x - 0.5, -4x + 3.5) ,0.0,1.0*br),
              clamp(min(4x + 0.5, -4x + 2.5) ,0.0,1.0*br))  
	      for x in LinRange(0.0,1.0,12)]

    mesh_col = mesh(points, jet[mesh_col])
    setmesh(mesh_col)
    poly(points, :fill, close=true)
    setline(1)
    sethue("black")
    poly(points, :stroke, close=true)
end

function bounds3D(mesh)
  max_x = []
  min_x = []
  max_y = []
  min_y = []
  max_z = []
  min_z = []
  for key in keys(mesh.nodes)
    if isempty(max_x)
      max_x = mesh.nodes[key][1]
      min_x = mesh.nodes[key][1]
      max_y = mesh.nodes[key][2]
      min_y = mesh.nodes[key][2]
      max_z = mesh.nodes[key][3]
      min_z = mesh.nodes[key][3]
    end
    if mesh.nodes[key][1] > max_x
      max_x = mesh.nodes[key][1]
    end    
    if mesh.nodes[key][1] < min_x
      min_x = mesh.nodes[key][1]
    end
    if mesh.nodes[key][1] > max_y
      max_y = mesh.nodes[key][2]
    end    
    if mesh.nodes[key][1] < min_y
      min_y = mesh.nodes[key][2]
    end
    if mesh.nodes[key][1] > max_z
      max_z = mesh.nodes[key][3]
    end    
    if mesh.nodes[key][1] < min_z
      min_z = mesh.nodes[key][3]
    end
  end 
  return max_x, min_x, max_y, min_y, max_z, min_z
end

function triangle_normal(P1,P2,P3)
  #Calculating normal from three points.
  x1,y1,z1 = P1
  x2,y2,z2 = P2
  x3,y3,z3 = P3
  n = [(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1)
       (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)
       (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)] 
  return n
end

function sorted_tri(mesh)
  #Counting all the Tet10 elements
  nTet10 = count(i->(i==:Tet10),values(mesh.element_types))
  #Creating an array for to keep all nodal data and normal for the side triangles
  triangles = zeros(nTet10*4,9+3)
  i = 1
  #Collecting the element data
  for key in keys(mesh.elements)
    if mesh.element_types[key] == :Tet10
      nodes = mesh.elements[key]
      xyz1 = mesh.nodes[nodes[1]]
      xyz2 = mesh.nodes[nodes[2]]
      xyz3 = mesh.nodes[nodes[3]]
      xyz4 = mesh.nodes[nodes[4]]
      #1st triangle
      nodes = [xyz1[1],xyz1[2],xyz1[3],xyz2[1],xyz2[2],xyz2[3],xyz4[1],xyz4[2],xyz4[3]]
      nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
      triangles[i,:] = vcat(nodes,nTri)
      i += 1
      #2nd triangle
      nodes = [xyz2[1],xyz2[2],xyz2[3],xyz3[1],xyz3[2],xyz3[3],xyz4[1],xyz4[2],xyz4[3]]
      nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
      triangles[i,:] = vcat(nodes,nTri)
      i += 1
      #3rd triangle
      nodes = [xyz3[1],xyz3[2],xyz3[3],xyz1[1],xyz1[2],xyz1[3],xyz4[1],xyz4[2],xyz4[3]]
      nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
      triangles[i,:] = vcat(nodes,nTri)
      i += 1
      #4th triangle
      nodes = [xyz1[1],xyz1[2],xyz1[3],xyz3[1],xyz3[2],xyz3[3],xyz2[1],xyz2[2],xyz2[3]]
      nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
      triangles[i,:] = vcat(nodes,nTri)
      i += 1
     end
  end
  return triangles
end

function trias_h5(h5data)
  #Counting all the Tet10 elements
  nTet10 = length(h5data["Topology"]["Tet10"]["Element IDs"])
  #Creating an array for to keep all nodal data and normal for the side triangles
  triangles = zeros(nTet10*4,9+3)
  i = 1
  #Collecting the element data
  for key in range(1,nTet10)
    nodes = h5data["Topology"]["Tet10"]["Connectivity"][:,key]

    xyz1 = h5data["Geometry"][:,nodes[1]+1]
    xyz2 = h5data["Geometry"][:,nodes[2]+1]
    xyz3 = h5data["Geometry"][:,nodes[3]+1]
    xyz4 = h5data["Geometry"][:,nodes[4]+1]

    #1st triangle
    nodes = [xyz1[1],xyz1[2],xyz1[3],xyz2[1],xyz2[2],xyz2[3],xyz4[1],xyz4[2],xyz4[3]]
    nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
    triangles[i,:] = vcat(nodes,nTri)
    i += 1
    #2nd triangle
    nodes = [xyz2[1],xyz2[2],xyz2[3],xyz3[1],xyz3[2],xyz3[3],xyz4[1],xyz4[2],xyz4[3]]
    nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
    triangles[i,:] = vcat(nodes,nTri)
    i += 1
    #3rd triangle
    nodes = [xyz3[1],xyz3[2],xyz3[3],xyz1[1],xyz1[2],xyz1[3],xyz4[1],xyz4[2],xyz4[3]]
    nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
    triangles[i,:] = vcat(nodes,nTri)
    i += 1
    #4th triangle
    nodes = [xyz1[1],xyz1[2],xyz1[3],xyz3[1],xyz3[2],xyz3[3],xyz2[1],xyz2[2],xyz2[3]]
    nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
    triangles[i,:] = vcat(nodes,nTri)
    i += 1

  end
  return triangles
end

function pproject(rP,r0,n,e1,e2)
  #Function projects 3D points to 2D plane with a dot product
  #Returns, 2D coordinates and distance from plane
  g = rP-r0
  s = n'g
  t1 = e1'g
  t2 = e2'g
  return [t1,t2,s]
end

function rangeInt(x,xMin,xMax,n)
  vn = collect(xMin:((xMax-xMin)/n):xMax)
  ff = searchsorted(vn,x)
  idx = ff.stop
  if idx > n
    idx = n
  end
  if idx < 1
    idx = 1
  end  
  return idx
end

function draw_triangles(triangles,plane,imageSize, deform, h5data)
  if plane=="XY"
    r0 = [0.,0.,0.]
    e1 = [1.,0.,0.]
    e2 = [0.,1.,0.]
    n  = [0.,0.,1.]
  end
  if plane=="YZ"
    r0 = [0.,0.,0.]  
    e1 = [0.,1.,0.]
    e2 = [0.,0.,1.]
    n  = [1.,0.,0.]
  end
  if plane=="XZ"
    r0 = [0.,0.,0.]  
    e1 = [1.,0.,0.]
    e2 = [0.,0.,1.]
    n  = [0.,1.,0.]
  end
  if plane=="XYZ"
    r0 = [0.,0.,0.]  
    e1 = [-1/sqrt(2),1/sqrt(2),0.]
    e2 = [0.,0.,1.]
    n  = [1/sqrt(3),1/sqrt(3),1/sqrt(3)]
  end
  
  if deform>-10000.0
    #Counting all the Tet10 elements
    nTet10 = length(h5data["Topology"]["Tet10"]["Element IDs"])
    #Creating an array for to keep all nodal data and normal for the side triangles
    triangles = zeros(nTet10*4,9+3)
    i = 1
    #Collecting the element data
    for key in range(1,nTet10)
      nodes = h5data["Topology"]["Tet10"]["Connectivity"][:,key]
      
      modeNodes=h5data["Results"]["Natural Frequency Analysis"]["Displacement"]["Mode 1"]
      
      xyz1 = h5data["Geometry"][:,nodes[1]+1] + modeNodes[:,nodes[1]+1]*deform 
      xyz2 = h5data["Geometry"][:,nodes[2]+1] + modeNodes[:,nodes[2]+1]*deform
      xyz3 = h5data["Geometry"][:,nodes[3]+1] + modeNodes[:,nodes[3]+1]*deform
      xyz4 = h5data["Geometry"][:,nodes[4]+1] + modeNodes[:,nodes[4]+1]*deform

      #1st triangle
      nodes = [xyz1[1],xyz1[2],xyz1[3],xyz2[1],xyz2[2],xyz2[3],xyz4[1],xyz4[2],xyz4[3]]
      nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
      triangles[i,:] = vcat(nodes,nTri)
      i += 1
      #2nd triangle
      nodes = [xyz2[1],xyz2[2],xyz2[3],xyz3[1],xyz3[2],xyz3[3],xyz4[1],xyz4[2],xyz4[3]]
      nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
      triangles[i,:] = vcat(nodes,nTri)
      i += 1
      #3rd triangle
      nodes = [xyz3[1],xyz3[2],xyz3[3],xyz1[1],xyz1[2],xyz1[3],xyz4[1],xyz4[2],xyz4[3]]
      nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
      triangles[i,:] = vcat(nodes,nTri)
      i += 1
      #4th triangle
      nodes = [xyz1[1],xyz1[2],xyz1[3],xyz3[1],xyz3[2],xyz3[3],xyz2[1],xyz2[2],xyz2[3]]
      nTri = triangle_normal(nodes[1:3],nodes[4:6],nodes[7:9])
      triangles[i,:] = vcat(nodes,nTri)
      i += 1
    end
  end
  
  #Projecting and calculating distance and angle from viewplane
  projTri = zeros(size(triangles,1),8)

  for  i=1:size(triangles,1)
    x1, y1, d1 = pproject(triangles[i,1:3],r0,n,e1,e2)
    x2, y2, d2 = pproject(triangles[i,4:6],r0,n,e1,e2)
    x3, y3, d3 = pproject(triangles[i,7:9],r0,n,e1,e2)
    d = (d1+d2+d3)/3
    nt = triangles[i,10:12]
    angle = acos((n'nt)/sqrt(nt'nt+n'n))
    projTri[i,:] = [x1,y1,x2,y2,x3,y3,d,angle]
  end
  
  #Centering xy-data
  xMax = maximum(projTri[:,[1,3,5]])
  xMin = minimum(projTri[:,[1,3,5]])
  yMax = maximum(projTri[:,[2,4,6]])
  yMin = minimum(projTri[:,[2,4,6]])
  
  projTri[:,1] = projTri[:,1].-(xMin+xMax)/2 
  projTri[:,2] = projTri[:,2].-(yMin+yMax)/2 
  projTri[:,3] = projTri[:,3].-(xMin+xMax)/2 
  projTri[:,4] = projTri[:,4].-(yMin+yMax)/2 
  projTri[:,5] = projTri[:,5].-(xMin+xMax)/2 
  projTri[:,6] = projTri[:,6].-(yMin+yMax)/2 
  
  #Scaling to figure size
  xMax = maximum(projTri[:,[1,3,5]])
  xMin = minimum(projTri[:,[1,3,5]])
  yMax = maximum(projTri[:,[2,4,6]])
  yMin = minimum(projTri[:,[2,4,6]])
  
  scaleX=imageSize[1]/(xMax-xMin)
  scaleY=imageSize[2]/(yMax-yMin)
  scale = minimum([scaleX,scaleY])*0.95
  
  projTri[:,1:6] = projTri[:,1:6]*scale
  
  projTri = sortslices(projTri, dims=1, by=x->x[7])
  
  #Color range index for plotting
  valueMax = maximum(projTri[:,[2,4,6]])
  valueMin = minimum(projTri[:,[2,4,6]])
  nValue = 12
  
  #Plotting the trianguladed surfaces with color scale
  for i=1:size(projTri,1)
    if projTri[i,8]<=(pi/2)
      nCol = [rangeInt(projTri[i,2],valueMin,valueMax,nValue),
	      rangeInt(projTri[i,4],valueMin,valueMax,nValue),
	      rangeInt(projTri[i,6],valueMin,valueMax,nValue)]
      br = (pi/2-projTri[i,8])/(pi/2)
      triangle([Point(projTri[i,1], projTri[i,2]),
		Point(projTri[i,3], projTri[i,4]),
		Point(projTri[i,5], projTri[i,6])],nCol,br)
    end
  end
end
