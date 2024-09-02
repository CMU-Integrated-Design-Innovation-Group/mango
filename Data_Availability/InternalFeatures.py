import trimesh

data = trimesh.load(r"C:\Users\Anthony\Documents\GitHub\Testing\input\test2.ply")
data.fix_normals()
data.export(r"C:\Users\Anthony\Documents\GitHub\Testing\input\test3.ply", encoding='ascii')