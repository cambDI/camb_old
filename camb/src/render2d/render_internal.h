/****************************************************************************
 * Copyright (C) 2009-2012 GGA Software Services LLC
 *
 * This file is part of Indigo toolkit.
 *
 * This file may be distributed and/or modified under the terms of the
 * GNU General Public License version 3 as published by the Free Software
 * Foundation and appearing in the file LICENSE.GPL included in the
 * packaging of this file.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 ***************************************************************************/

#ifndef __render_internal_h__
#define __render_internal_h__

#include "render_common.h"

namespace indigo {

class RenderContext;

class MoleculeRenderInternal {
public:
   MoleculeRenderInternal (const RenderOptions& opt, const RenderSettings& settings, RenderContext& cw);
   void setMolecule (BaseMolecule* mol);
   void setIsRFragment (bool isRFragment);
   void setScaleFactor (const float scaleFactor, const Vec2f& min, const Vec2f& max);
   void render ();

   void setReactionComponentProperties (const Array<int>* aam, const Array<int>* reactingCenters, const Array<int>* inversions);
   void setQueryReactionComponentProperties (const Array<int>* exactChanges);

   DECL_ERROR;
private:
   enum STEREOGROUPS_MODE {STEREOGROUPS_SHOW, STEREOGROUPS_HIDE};
   struct LocalOptions {
      STEREOGROUPS_MODE stereoMode;
   };

   BondEnd& _be (int beid);
   const BondEnd& _be (int beid) const;
   BondDescr& _bd (int bid);
   const BondDescr& _bd (int bid) const;
   AtomDesc& _ad (int aid);
   const AtomDesc& _ad (int aid) const;
   void _checkSettings ();
   void _extendRenderItem (RenderItem& item, const float extent);
   bool _clipRaySegment (float& offset, const Vec2f& p, const Vec2f& d, const Vec2f& n0, const Vec2f& a, const Vec2f& b, const float w);
   bool _clipRayBox (float& offset, const Vec2f& p, const Vec2f& d, const Vec2f& rp, const Vec2f& sz, const float w);
   void _findMinMax();
   void _objCoordTransform(Vec2f& p, const Vec2f& v) const;
   void _objDistTransform(Vec2f& p, const Vec2f& v) const;
   void _initCoordinates();
   void _determineDoubleBondShift();
   void _determineStereoGroupsMode();
   static const char* _getStereoGroupText (int type);
   bool _ringHasSelfIntersectionsSimple(const Ring& ring);
   bool _ringHasSelfIntersections(const Ring& ring);
   void _findRings();
   void _prepareLabels();
   void _rotateHalfCenteredBonds();
   bool _isSingleHighlighted (int aid);
   bool _vertexIsHighlighted (int aid);
   bool _edgeIsHighlighted (int bid);
   bool _hasQueryModifiers (int aid);
   void _findNearbyAtoms();
   void _initHydroPos(int aid);
   int  _hydroPosFindConflict(int i);
   bool _hydroPosCorrectGreedy ();
   void _hydroPosCorrectRepulse ();
   void _initAtomData();
   void _initRGroups();
   void _loadBrackets(SGroup& sg, const Array<Vec2f[2]>& coord, bool transformCoordinates);
   void _placeBrackets(SGroup& sg, const Array<int>& atoms);
   void _positionIndex(SGroup& sg, int ti, bool lower);
   void _loadBracketsAuto(const BaseMolecule::SGroup& group, SGroup& sg);
   void _initDataSGroups();
   void _initSruGroups();
   void _initMulGroups();
   void _initSupGroups();
   void _prepareSGroups();
   void _findAnglesOverPi();
   void _renderBondIds();
   void _renderAtomIds();
   void _renderEmptyRFragment();
   void _renderLabels();
   void _renderRings();
   void _renderSGroups ();
   void _setHighlightOpt();
   void _resetHighlightOpt();
   void _renderBonds();
   void _applyBondOffset();
   void _setBondCenter ();
   float _getBondOffset (int aid, const Vec2f& pos, const Vec2f& dir, const float bondWidth);
   void _calculateBondOffset();
   void _findNeighbors();
   void _findCenteredCase();
   void _initBondData();
   void _initBondEndData();
   void _extendRenderItems();
   BondEnd& _getBondEnd(int aid, int nei);
   int _getBondEndIdx (int aid, int nei);
   int _getOpposite (int beid) const;
   void _drawAtom (const AtomDesc& desc);
   void _writeQueryAtomToString (Output& output, int aid);
   bool _writeDelimiter (bool needDelimiter, Output &output);
   void _writeQueryModifier (Output& output, int aid);
   int _findClosestCircle (Vec2f& p, int aid, float radius, int skip = -1);
   int _findClosestBox (Vec2f& p, int aid, const Vec2f& sz, float mrg, int skip = -1);
   void _preparePseudoAtom (int aid, int color, bool highlighted);
   void _prepareLabelText (int aid);
   void _prepareAAM ();
   int _pushTextItem (RenderItem::TYPE type, int color, bool highlighted);
   int _pushTextItem (AtomDesc& ad, RenderItem::TYPE type, int color, bool highlighted);
   int _pushTextItem (SGroup& sg, RenderItem::TYPE ritype, int color = CWC_BASE);
   int _pushGraphItem (RenderItem::TYPE type, int color, bool highlighted);
   int _pushGraphItem (AtomDesc& ad, RenderItem::TYPE type, int color, bool highlighted);
   int _pushGraphItem (SGroup& ad, RenderItem::TYPE type, int color = CWC_BASE);
   const char* _valenceText (const int valence);
   float _ctghalf (float cs);
   void _drawBond (int b);
   void _drawTopology (BondDescr& bd);
   void _drawReactingCenter (BondDescr& bd, int rc);
   float _doubleBondShiftValue (const BondEnd& be, bool right, bool centered);
   void _prepareDoubleBondCoords (Vec2f* coord, BondDescr& bd, const BondEnd& be1, const BondEnd& be2, bool allowCentered);
   void _drawStereoCareBox (BondDescr& bd, const BondEnd& be1, const BondEnd& be2);
   void _bondSingle (BondDescr& bd, const BondEnd& be1, const BondEnd& be2);
   void _bondDouble (BondDescr& bd, const BondEnd& be1, const BondEnd& be2);
   void _bondSingleOrAromatic (BondDescr& bd, const BondEnd& be1, const BondEnd& be2);
   void _bondDoubleOrAromatic (BondDescr& bd, const BondEnd& be1, const BondEnd& be2);
   void _bondSingleOrDouble (BondDescr& bd, const BondEnd& be1, const BondEnd& be2);
   void _bondAromatic (BondDescr& bd, const BondEnd& be1, const BondEnd& be2);
   void _bondTriple (BondDescr& bd, const BondEnd& be1, const BondEnd& be2);
   void _bondAny (BondDescr& bd, const BondEnd& be1, const BondEnd& be2);

   // local
   void* _hdc;
   BaseMolecule* _mol;
   RenderContext& _cw;
   float _scale;
   Vec2f _min, _max;
   LocalOptions _lopt;
   bool isRFragment;
   const RenderSettings& _settings;
   const RenderOptions& _opt;
   TL_CP_DECL(MoleculeRenderData, _data);
   TL_CP_DECL(Array<int>, _atomMapping);
   TL_CP_DECL(Array<int>, _atomMappingInv);
   TL_CP_DECL(BaseMolecule::Mapping, _bondMappingInv);
};

}

#endif //__render_internal_h__
