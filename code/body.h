#include"constants.h"
#include"vector.h"
#include"collider.h"

class Body{
    private:
        Vector3D r, V, F; double m, R;
    public:
        Body();
        friend class Collider;
};

Body::Body(){
    
}

