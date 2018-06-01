import bpy
from bpy.types import Operator
from bpy.types import Panel
from bpy.props import *

bl_info = {'name': 'Elfin Front', 'category': 'Object'}

# Addon design notes
#   * Each separate object is a separate chain
#   * There should be no faces in any object
#       * Code simply ignores faces
#       x Provide a face delete operator 
#   * There should be no discontinuities in an object
#       * Code should verify this
#   x Unit conversion: 1 blender unit is 10 A, or 1 nm
#       x In other words 1 unit in blender is 10 pymol units
#   * Module avatar generation:
#       * Singles are pairs exported as OBJ from Pymol
#           * Pre-process OBJs:
#               * Load OBJ
#               * Remove doubles (edit mode)
#                   bpy.ops.mesh.remove_doubles()
#               * Decimate to 0.05 ratio
#                   bpy.ops.object.modifier_add(type='DECIMATE')
#                   bpy.context.object.modifiers['Decimate'].ratio = 0.05
#                   bpy.ops.object.modifier_apply(apply_as='DATA', modifier='Decimate')
#               * Save as new OBJ
#       * How to do avatars in viewport elegantly?
#
# bpy.ops.elfin.unload()
# filename = 'C:\\Users\\Akaoni\\Desktop\\ElfinWork\\elfin\\scripts\\Blender\\ElfinFrontBlenderAddon.py'
# exec(compile(open(filename).read(), filename, 'exec'))

class ElfinPropertyGroup(bpy.types.PropertyGroup):
    PpSrcDir = StringProperty(subtype='DIR_PATH')
    PpDstDir = StringProperty(subtype='DIR_PATH')

class ElfinDebugPanel(Panel):
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_label = 'Debug'
    bl_context = 'objectmode'
    bl_category = 'Elfin'


    def draw(self, context):
        layout = self.layout
        row = layout.row(align=True)
        col = row.column()
        col.operator('elfin.unload', text='Unload')
        col.operator('elfin.delete_faces', text='Delete Faces')

class ElfinPreprocessPanel(Panel):
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_label = 'Preprocess'
    bl_context = 'objectmode'
    bl_category = 'Elfin'


    def draw(self, context):
        scn = context.scene

        layout = self.layout

        row = layout.row(align=True)
        col = row.column()
        col.prop(scn.ElfinFront, 'PpSrcDir', text='Source')
        col.prop(scn.ElfinFront, 'PpDstDir', text='Destination')
        col.operator('elfin.preproess_objs', text='Preprocess')

class ElfinExportPanel(Panel):
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

class ElfinPreprocessObjsOperator(Operator):
    bl_idname = 'elfin.preproess_objs'
    bl_label = 'Preprocess module obj files'

    
    def execute(self, context):
        print('Dialog?')
        return {'FINISHED'}
 
    def invoke(self, context, event):
        self.srcPath = ''
        return context.window_manager.invoke_props_dialog(self)

class ElfinExportOperator(Operator):
    bl_idname = 'elfin.export'
    bl_label = 'Export as Elfin input'

    def execute(self, context):
        # Each separate object is a separate chain
        print('Unimplemented')

        return {'FINISHED'}

class ElfinUnloadOperator(Operator):
    bl_idname = 'elfin.unload'
    bl_label = 'Unload Elfin Front'

    def execute(self, context):
        unregister()
        return {'FINISHED'}

class ElfinDeleteFacesOperator(Operator):
    bl_idname = 'elfin.delete_faces'
    bl_label = 'Delete Faces (selected only)'
    
    def execute(self, context):
        bpy.ops.object.mode_set(mode='OBJECT')
        for obj in bpy.context.selected_objects:
            bpy.context.scene.objects.active = obj
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.delete(type='ONLY_FACE')
            bpy.ops.object.mode_set(mode='OBJECT')
        return {'FINISHED'}

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