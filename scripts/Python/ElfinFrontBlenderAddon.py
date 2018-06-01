import bpy
from bpy.types import Operator

bl_info = {'name': 'Elfin Front Addon', 'category': 'Object'}

# Addon design notes
#   * Each separate object is a separate chain
#   * There should be no faces in any object
#       * Code simply ignores faces
#       x Provide a face delete operator 
#   * There should be no discontinuities in an object
#       * Code should verify this
#   x Unit conversion: 1 blender unit is 10 A, or 1 nm
#   * Module avatar generation:
#       * Singles are pairs exported as OBJ from Pymol
#           * Pre-process OBJs:
#               * Load OBJ
#               * Decimate to 0.05 ratio
#               * Save as new OBJ
#       * How to do avatars in viewport elegantly?

# API notes
# filename = 'C:\\Users\\Akaoni\\Desktop\\ElfinWork\\elfin\\scripts\\Python\\ElfinFrontBlenderAddon.py'; exec(compile(open(filename).read(), filename, 'exec'))
# bpy.ops.object.mode_set(mode='OBJECT')
# bpy.ops.object.select_all(action='SELECT')
# bpy.ops.object.editmode_toggle()
# bpy.ops.mesh.select_all(action='SELECT')
# bpy.ops.mesh.select_all(action='DESELECT')
# bpy.ops.mesh.delete(type='ONLY_FACE')



class ElfinExportOperator(Operator):
    bl_idname = 'elfin.elfin_export_operator'
    bl_label = 'Export as Elfin input'

    def execute(self, context):
        # Each separate object is a separate chain

        return {'FINISHED'}

class FaceDelOperator(Operator):
    bl_idname = 'elfin.face_delete_operator'
    bl_label = 'Delete Faces (selected only)'
    
    def execute(self, context):
        bpy.ops.object.mode_set(mode='OBJECT')
        for obj in bpy.context.selected_objects:
            bpy.context.scene.objects.active = obj
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.delete(type='ONLY_FACE')
            bpy.ops.object.mode_set(mode='OBJECT')
        return {'FINISHED'}

classesToRegister = [FaceDelOperator]

def register():
    for cls in classesToRegister:
        bpy.utils.register_class(cls)
    print('Elfin Front Addon registered')
    
def unregister():
    for cls in classesToRegister:
        bpy.utils.unregister_class(cls)
    print('Elfin Front Addon unregistered')
    
if __name__ == '__main__':
    register()