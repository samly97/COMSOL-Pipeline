HOME = pwd;
COMSOL_DIR = '/Applications/COMSOL56/Multiphysics/mli';

cd(COMSOL_DIR);
mphstart(2036);
cd(HOME);

import com.comsol.model.*
import com.comsol.model.util.*;

disp('LiveLink is ready')