classdef objectHandle < handle
   properties
      obj
   end
   methods
      function h_obj = objectHandle(obj)
          h_obj.obj = obj;
      end
   end
end

