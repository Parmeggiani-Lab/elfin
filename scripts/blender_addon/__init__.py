
#
# Elfin's Front end as a Blender addon
# 
# Author: Joy Yeh 
# Email: joyyeh.tw@gmail.com
#

import bpy
from bpy.props import *

import glob, os

bl_info = {'name': 'Elfin Front', 'category': 'Elfin'}

# Addon design notes
#   * Each separate object is a separate chain
#   * There should be no faces in any object
#       * Code simply ignores faces
#       x Provide a face delete operator 
#   * There should be no discontinuities in an object
#       * Code should verify this
#   x Unit conversion: 1 blender unit is 10 A, or 1 nm
#       x 1 blender unit === 10 pymol units
#   * Module avatar generation:
#       * Singles and pairs exported as OBJ from Pymol
#           * Pre-process OBJs
#       * How to do avatars in viewport elegantly?
#

class ElfinPropertyGroup(bpy.types.PropertyGroup):
  pp_src_dir = bpy.props.StringProperty(subtype='DIR_PATH', default='C:\\Users\\Akaoni\\Desktop\\ElfinWork\\elfin\\resources\\obj_aligned\\')
  pp_dst_dir = bpy.props.StringProperty(subtype='FILE_PATH', default='C:\\Users\\Akaoni\\Desktop\\ElfinWork\\elfin\\resources\\elfin_front_library.blend')
  pp_decimate_ratio = bpy.props.FloatProperty(default=0.15, min=0.00, max=1.00)


# -------------------- Panels --------------------

class ElfinDebugPanel(bpy.types.Panel):
  bl_space_type = 'VIEW_3D'
  bl_region_type = 'TOOLS'
  bl_label = 'Debug'
  bl_context = 'objectmode'
  bl_category = 'Elfin'

  def draw(self, context):
    layout = self.layout
    row = layout.row(align=True)
    col = row.column()
    col.operator('elfin.reset', text='Reset Properties')
    col.operator('elfin.delete_faces', text='Delete faces (selection)')
    col.operator('elfin.load_all_obj_files', text='Load all objects')
    col.operator('elfin.process_obj', text='Process object (selection)')

class ElfinProcessPanel(bpy.types.Panel):
  bl_space_type = 'VIEW_3D'
  bl_region_type = 'TOOLS'
  bl_label = 'Process'
  bl_context = 'objectmode'
  bl_category = 'Elfin'

  def draw(self, context):
    scn = context.scene
    layout = self.layout

    row = layout.row(align=True)
    col = row.column()
    col.prop(scn.ElfinFront, 'pp_src_dir', text='Source')
    col.prop(scn.ElfinFront, 'pp_dst_dir', text='Destination')
    col.prop(scn.ElfinFront, 'pp_decimate_ratio', text='Decimate Ratio')
    col.operator('elfin.process_all', text='Process')

class ElfinExportPanel(bpy.types.Panel):
  bl_space_type = 'VIEW_3D'
  bl_region_type = 'TOOLS'
  bl_label = 'Export'
  bl_context = 'objectmode'
  bl_category = 'Elfin'

  def draw(self, context):
    layout = self.layout
    row = layout.row(align=True)
    col = row.column()
    col.operator('elfin.export', text='Export to .ei')


# -------------------- Operators --------------------

class ElfinExportOperator(bpy.types.Operator):
  bl_idname = 'elfin.export'
  bl_label = 'Export as Elfin input'

  def execute(self, context):
    # Each separate object is a separate chain
    print('Unimplemented')

    return {'FINISHED'}

def make_dir(path):
  try:
    os.makedirs(path)
  except OSError as exc:  # Python >2.5
    if os.path.isdir(path):
      pass
    else:
      raise

class ElfinLoadAllObjFiles(bpy.types.Operator):
  bl_idname = 'elfin.load_all_obj_files'
  bl_label = 'Load all Pymol obj files'

  def execute(self, context):
    abs_src_path = bpy.path.abspath(bpy.context.scene.ElfinFront.pp_src_dir)
    src_path_list = glob.glob(abs_src_path + '*')
    src_folders = list(map(os.path.basename, src_path_list))

    module_types = ['singles', 'doubles', 'hubs']
    
    if sum([folder in src_folders for folder in module_types]) != 3:
      self.report({'ERROR'}, 'Source folder {} does not contain {} folders: {}'. \
        format(abs_src_path, module_types, src_folders))
      return {'CANCELLED'}

    obj_files = [f for flist in [glob.glob(abs_src_path + '/{}/*.obj'.format(mt)) for mt in module_types] for f in flist]

    for src_obj_file in obj_files:
      bpy.ops.import_scene.obj(filepath=src_obj_file)
    
    return {'FINISHED'}

class ElfinProcessAllOperator(bpy.types.Operator):
  bl_idname = 'elfin.process_all'
  bl_label = 'Process module obj files'

  def execute(self, context):
    if len(bpy.context.scene.objects) != 0:
      self.report({'ERROR'}, 'Scene must be empty for processing')
      return {'CANCELLED'}
    
    bpy.ops.elfin.load_all_obj_files()

    make_dir(bpy.path.dirname(bpy.context.scene.ElfinFront.pp_dst_dir))
    abs_dst_path = bpy.path.abspath(bpy.context.scene.ElfinFront.pp_dst_dir)

    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.elfin.process_obj()
    bpy.ops.wm.save_as_mainfile(filepath=abs_dst_path, copy=True)
    
    return {'FINISHED'}

class ElfinProcessObjOperator(bpy.types.Operator):
  bl_idname = 'elfin.process_obj'
  bl_label = 'Process module object'

  def execute(self, context):
    for obj in bpy.context.selected_objects:
      # Shrink to scale
      obj.scale = (.1, .1, .1)
      bpy.context.scene.objects.active = obj
      
      # Fix normals and remove superimposed vertices. Do this before
      # decimate so the ratio works as intended.
      bpy.ops.object.mode_set(mode='EDIT')
      bpy.ops.mesh.customdata_custom_splitnormals_clear()
      bpy.ops.mesh.remove_doubles()

      # Reduce polygons
      bpy.ops.object.mode_set(mode='OBJECT')
      bpy.ops.object.modifier_add(type='DECIMATE')
      bpy.context.object.modifiers['Decimate'].ratio = context.scene.ElfinFront.pp_decimate_ratio
      bpy.ops.object.modifier_apply(apply_as='DATA', modifier='Decimate')
    return {'FINISHED'}

class ElfinDeleteFacesOperator(bpy.types.Operator):
  bl_idname = 'elfin.delete_faces'
  bl_label = 'Delete Faces (selected only)'
  
  def execute(self, context):
    selObjs = bpy.context.selected_objects
    for obj in selObjs:
      bpy.context.scene.objects.active = obj
      bpy.ops.object.mode_set(mode='EDIT')
      bpy.ops.mesh.delete(type='ONLY_FACE')
      bpy.ops.object.mode_set(mode='OBJECT')
    return {'FINISHED'}

class ElfinResetOperator(bpy.types.Operator):
  bl_idname = 'elfin.reset'
  bl_label = 'Reset Elfin Front properties'

  def execute(self, context):
    scn = context.scene
    for k in scn.ElfinFront.keys():
      scn.ElfinFront.property_unset(k)
    return {'FINISHED'}

# -------------------- Addon Basic Functions --------------------

def register():
  bpy.utils.register_module(__name__)
  bpy.types.Scene.ElfinFront = bpy.props.PointerProperty(type=ElfinPropertyGroup)
  print('Elfin Front Addon registered')
  
def unregister():
  bpy.utils.unregister_module(__name__)
  del bpy.types.Scene.ElfinFront
  print('Elfin Front Addon unregistered')
  
if __name__ == '__main__':
  register()