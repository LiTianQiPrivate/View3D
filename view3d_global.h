
#ifndef VIEW3D_GLOBAL_H
#define VIEW3D_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(VIEW3D_LIBRARY)
#  define VIEW3D_EXPORT Q_DECL_EXPORT
#else
#  define VIEW3D_EXPORT Q_DECL_IMPORT
#endif

#endif // VIEW3D_GLOBAL_H
