// SPDX-License-Identifier: GPL-3.0-or-later
#include "thicknessPolarFamily.H"
#include "primitiveEntry.H"

#include <algorithm>
#include <cmath>
#include <utility>

Foam::thicknessPolarFamily::thicknessPolarFamily
(
    const dictionary& blade,
    const dictionary& profiles,
    const wordList& names
)
:
    profiles_(profiles),
    names_(names)
{
    if (!blade.readIfPresent("profileThicknesses", thicknesses_)) return;

    if (thicknesses_.empty() || thicknesses_.size() != names_.size())
    {
        FatalIOErrorInFunction(blade)
            << "profileThicknesses requires one t/c fraction per elementProfiles entry"
            << exit(FatalIOError);
    }

    tables_.setSize(names_.size());
    forAll(thicknesses_, i)
    {
        const scalar tc = thicknesses_[i];
        if (!(tc > 0 && tc <= 1) || (i && tc <= thicknesses_[i-1]))
        {
            FatalIOErrorInFunction(blade)
                << "profileThicknesses must be strictly increasing fractions in (0, 1]"
                << exit(FatalIOError);
        }
        const dictionary& profile = profiles_.subDict(names_[i]);
        if
        (
            profile.lookupOrDefault<word>("tableType", "singleRe") != "singleRe"
         || profile.lookupOrDefault<scalar>("Re", 0) > VSMALL
        )
        {
            FatalIOErrorInFunction(profile)
                << "Thickness blending requires static singleRe profiles without Re correction"
                << exit(FatalIOError);
        }

        List<scalarList>& table = tables_[i];
        profile.lookup("data") >> table;
        if (table.size() < 2 || table.size() != tables_[0].size())
        {
            FatalIOErrorInFunction(profile)
                << "The static polars must share an angle grid with at least two rows"
                << exit(FatalIOError);
        }
        forAll(table, j)
        {
            const scalarList& row = table[j];
            if
            (
                row.size() != 4
             || !std::all_of(row.begin(), row.end(), [](scalar v){ return std::isfinite(v); })
             || (j && row[0] <= table[j-1][0])
             || (i && row[0] != tables_[0][j][0])
            )
            {
                FatalIOErrorInFunction(profile)
                    << "Expected finite (alpha Cl Cd Cm) rows on one strictly increasing angle grid"
                    << exit(FatalIOError);
            }
        }
    }
}


Foam::dictionary Foam::thicknessPolarFamily::at
(
    const scalar thickness,
    word& profileName
) const
{
    if (!active() || !std::isfinite(thickness))
    {
        FatalIOErrorInFunction(profiles_)
            << "A thickness lookup requires a family and a finite element t/c"
            << exit(FatalIOError);
    }
    const label n = thicknesses_.size();
    const label upper = std::lower_bound(thicknesses_.begin(), thicknesses_.end(), thickness)
                      - thicknesses_.begin();
    if (upper == 0 || upper == n || thickness == thicknesses_[upper])
    {
        profileName = names_[min(upper, n-1)];
        return profiles_.subDict(profileName);
    }

    const label lower = upper - 1;
    const scalar w = (thickness - thicknesses_[lower])
                   / (thicknesses_[upper] - thicknesses_[lower]);
    // Passing a numeric list to dictionary::add formats it at the default
    // stream precision. Store scalar tokens directly to retain all digits.
    List<token> data(2 + 6*tables_[lower].size());
    label pos = 0;
    data[pos++] = token::BEGIN_LIST;
    forAll(tables_[lower], row)
    {
        data[pos++] = token::BEGIN_LIST;
        data[pos++] = tables_[lower][row][0];
        for (label col = 1; col <= 3; ++col)
        {
            data[pos++] = (1-w)*tables_[lower][row][col] + w*tables_[upper][row][col];
        }
        data[pos++] = token::END_LIST;
    }
    data[pos] = token::END_LIST;
    profileName = names_[lower] + "_" + names_[upper];
    dictionary result;
    result.add("tableType", word("singleRe"));
    result.add(new primitiveEntry("data", std::move(data)));
    return result;
}
