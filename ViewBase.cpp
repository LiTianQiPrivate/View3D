#include "ViewBase.h"

ViewBase::ViewBase()
{

}

ViewBase::MOUSE_ACTION ViewBase::getMouseAction() const
{
    return mouseAction;
}

void ViewBase::setMouseAction(const MOUSE_ACTION &value)
{
    mouseAction = value;
}
