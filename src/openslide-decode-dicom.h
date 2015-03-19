/*
 *  OpenSlide, a library for reading whole slide image files
 *
 *  Copyright (c) 2015 Mathieu Malaterre
 *  All rights reserved.
 *
 *  OpenSlide is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, version 2.1.
 *
 *  OpenSlide is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with OpenSlide. If not, see
 *  <http://www.gnu.org/licenses/>.
 *
 */

#ifndef OPENSLIDE_OPENSLIDE_DECODE_DICOM_H_
#define OPENSLIDE_OPENSLIDE_DECODE_DICOM_H_

#include "openslide-private.h"
#include "openslide-hash.h"

#include <stdio.h>
#include <stdint.h>
#include <glib.h>

/* DICOM container support */
/* Thread-safe. */

struct _openslide_dicom *_openslide_dicom_create(const char *filename,
                                                       GError **err);

void _openslide_dicom_destroy(struct _openslide_dicom *d);

bool _openslide_dicom_readindex(struct _openslide_dicom *instance, GError **err);

// helpful printout?
//void _openslide_dicomdir_print(struct _openslide_dicom *d);

#endif
