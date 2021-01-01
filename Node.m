classdef Node < matlab.mixin.Copyable
    
   properties
      id
      % Row id of original dataset that traveled here
      sample_ids
      decision_feature
      decision_val
      
      node_mean
      node_var
      
      linear_reg_coef
      ancestral_decision_features 
      % Ancestral nodes decision_features if this is a leaf node.
      % Correponding position in linear_reg_coef
   end
   
   methods
       % Initialize
       function obj = Node(id, decision_feature, decision_val)
           if nargin == 1
               obj.id = id;
               obj.decision_feature = NaN;
               obj.decision_val = NaN;
           elseif nargin == 3
               obj.id = id;
               obj.decision_feature = decision_feature;
               obj.decision_val = decision_val;
           end
           
           
               
       end
               
           
    end
       
       
end
    
