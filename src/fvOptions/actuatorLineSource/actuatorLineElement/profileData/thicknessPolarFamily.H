// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef thicknessPolarFamily_H
#define thicknessPolarFamily_H

#include "dictionary.H"
#include "scalarList.H"
#include "wordList.H"

namespace Foam
{

// Static section polars on a common angle grid, ordered by increasing t/c.
class thicknessPolarFamily
{
    const dictionary& profiles_;
    const wordList& names_;
    scalarList thicknesses_;
    List<List<scalarList>> tables_;

public:
    thicknessPolarFamily
    (
        const dictionary& blade,
        const dictionary& profiles,
        const wordList& names
    );

    bool active() const { return !thicknesses_.empty(); }

    // Build the existing singleRe dictionary once, at element creation.
    dictionary at(const scalar thickness, word& profileName) const;
};

} // End namespace Foam

#endif
