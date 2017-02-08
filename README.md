# contact_models_and_data
# data formating
The data in data/ellipse uniform is formated such that the first 6 elements of bounce_array(trial).states
are the configurations (x,y,theta) and (dx,dy,dtheta) pre-impact and the next 6 are post-impact measured,
all related to the center of mass of the object.

bounce_array(trial).n/.d are the normal and tangent Jacobians, applying these vectors will transform the velocity
of the COM to the contact point.
