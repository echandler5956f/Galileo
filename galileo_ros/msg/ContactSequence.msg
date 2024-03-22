# It is assumed that the contact surfaces are set in the order of the end effectors. 
# contact_surface_ids[i][j] = contact surface at phase i, and for end effector j. 
# If contact_surface_ids[i][j] = NO_SURFACE, EE is not in contact.
int32[][] contact_surface_ids

int32[] knot_nums
float64[] knot_times