"""
Module with utility data and functions.
"""

# Import Python standard libraries
import hashlib
import itertools
import random
from typing import Hashable, List, Optional

# Import 3rd-party libraries
import numpy as np


# Define data for random label generation: sound classes, too complex
# clusters, and syllable patterns
_SOUNDS = {
    "C": "b p d t f v s z r l g k m n h".split(),
    "V": "a e i o u".split(),
}
_COMPLEX_CLUSTERS = [
    "pb",
    "bp",
    "sz",
    "zs",
    "dl",
    "gk",
    "kg",
    "bd",
    "db",
    "zp",
    "pv",
    "pf",
    "sr",
]
_PATTERNS = ["V", "CV", "CV", "CVC"]


def set_seeds(seed: Optional[Hashable] = None):
    """
    Set seeds globally from the user provided one.

    The function takes care of reproducibility and allows to use strings and
    floats as seed for `numpy` as well.
    """

    # Set Python generator seed
    random.seed(seed)

    # Allows using strings as np seeds, which only takes uint32 or arrays of
    # numbers. As this does not allow to accept `None`, we manually
    # call the function if necessary,
    if seed is None:
        np.random.seed()
    else:
        if isinstance(seed, (str, float)):
            np_seed = np.frombuffer(
                hashlib.sha256(str(seed).encode("utf-8")).digest(), dtype=np.uint32
            )
        else:
            np_seed = seed

        # Set the np seed
        np.random.seed(np_seed)


def _gen_syl(min_syl: int, max_syl: int) -> List[str]:
    """
    Internal function for generating a random syllable.

    :param min_syl: The minimum number of syllables.
    :param max_syl: The maximum number of syllables.
    :return: A list of strings, each with one random syllable.
    """

    # Map each syllable to random sounds
    syllables = []
    for _ in range(np.random.randint(min_syl, max_syl + 1)):
        syllable = "".join(
            [
                np.random.choice(_SOUNDS[sound_class])
                for sound_class in np.random.choice(_PATTERNS)
            ]
        )

        syllables.append(syllable)

    return syllables


def _clean_label(label: str) -> str:
    """
    Returns a cleaned version of a label.

    This basically takes care of generating a more readable random
    label, removing harder to read geminates, treating "h" differently,
    etc. It is also used to guarantee that all labels will be
    capitalized.

    :param label: The label to be cleaned.
    :return: The cleaned version of the label.
    """

    # We "uncapitalize" label, as it might have been capitalized before
    # and having everything in lowercase makes the code easier to follow.
    label = label.lower()

    # Remove all "h" next to another consonant (including "h" itself),
    # making sure we only have "h" in intervocalic position or at the
    # beginning of the word (increases readability)
    label = label.replace("hh", "")
    for cons in _SOUNDS["C"]:
        label = label.replace(cons + "h", cons)
        label = label.replace("h" + cons, cons)

    # Remove clusters that are too complex by selecting one random sound
    for cluster in _COMPLEX_CLUSTERS:
        if cluster in label:
            label = label.replace(cluster, cluster[np.random.randint(0, 1)])

    # Remove geminated vowels
    for vowel in _SOUNDS["V"]:
        label = label.replace(vowel + vowel, vowel)

    # Replace initial "i" with "wi" -- this makes reading labels easier
    # in most typefaces.
    if label.startswith("i"):
        label = "w" + label

    return label.capitalize()


def random_labels(size: int = 1, seed: Optional[Hashable] = None) -> List[str]:
    """
    Returns a list of unique random pronounceable labels.

    :param size: The number of labels in the returned set. Defaults to 1.
    :param seed: An optional seed for the random number generator. Defaults
        to `None`.
    :return: The list of unique labels.
    """

    # Initialize the RNG
    set_seeds(seed)

    # Iterate until enough unique labels have been collected
    ret_labels = []
    for _ in range(size):
        # Generate a random capitalized label with 2 to 3 syllables.
        syllables = _gen_syl(2, 3)
        label = _clean_label("".join(syllables))

        # Append more syllables if necessary, one at a time, until an unique
        # name is generated.
        while True:
            if label not in ret_labels:
                break

            label = _clean_label(label + _gen_syl(1, 1)[0])

        # Collect the generated label.
        ret_labels.append(label)

    # Return the list of labels
    return ret_labels


# Please note that this library was originally relying on the more
# complex word-generator implemented in Abzu, with better results and better
# code. Still, it was fun to take a break from serious work and work in
# this pseudo-modern Latin, and we only really care about good-enough
# labels in this case. In other words, please don't care too much about
# this function or the quality of its code ;)
def random_species(size: int = 1, seed: Optional[Hashable] = None) -> List[str]:
    """
    Returns a list of unique random species labels.

    :param size: The number of labels in the returned set. Defaults to 1.
    :param seed: An optional seed for the random number generator. Defaults
        to `None`.
    :return: The list of unique labels.
    """

    # Obtain random labels: we need double the number of items, as
    # the first half of the list will be the genera,
    # the second the epithets (they will be combined later, with
    # capitalization, etc.)
    labels = [label.lower() for label in random_labels(size * 2, seed)]

    # Remove all the "h" in the original labels, as they are intended to be
    # IPA /h/ and not aspiration
    labels = [label.replace("h", "") for label in labels]

    # Replace all "f" with "ph" and "k" with "c" ("latinizing")
    labels = [label.replace("f", "ph") for label in labels]
    labels = [label.replace("k", "c") for label in labels]

    # If the label does not end in a vowel or in s/r, add a random suffix
    labels = [
        label + np.random.choice(_SOUNDS["V"]) + np.random.choice(["s", ""])
        if label[-1] not in _SOUNDS["V"] + ["s", "r"]
        else label
        for label in labels
    ]

    # All the labels ending in "u/e" will end in "us/es"
    labels = [
        label + "s" if (label.endswith("u") or label.endswith("e")) else label
        for label in labels
    ]

    # If a label ends in "i", it will end in "is" 75% of the time, otherwise
    # in "ii"
    labels = [
        label + "s" if (label.endswith("i") and np.random.random() <= 0.75) else label
        for label in labels
    ]
    labels = [label + "i" if label.endswith("i") else label for label in labels]

    # If a label ends in "a", it will end in "as" 50% of the time
    labels = [
        label + "s" if (label.endswith("a") and np.random.random() <= 0.5) else label
        for label in labels
    ]

    # In case of t/p in front of a vowel,
    # there is a 50% of change that it will be aspirated
    # (i.e., plus an "h"). As all the "h" at this point are intervocalic,
    # we can just run some replacements.
    for plosive in "tp":
        for vowel in _SOUNDS["V"]:
            labels = [
                label.replace(plosive + vowel, plosive + "h" + vowel)
                if np.random.random() <= 0.5
                else label
                for label in labels
            ]

    # In case the label starts with a labial plosive, there is 50% of chance
    # that it will gain an "s" in front
    labels = [
        "s" + label if (label.startswith("p") or label.startswith("b")) else label
        for label in labels
    ]

    # Increase the number of geminates in case of intervocalic consonants
    for cons in "bpdtsrlgmn":
        for vowel1 in _SOUNDS["V"]:
            for vowel2 in _SOUNDS["V"]:
                source = vowel1 + cons + vowel2
                target = vowel1 + cons + cons + vowel2

                labels = [
                    label.replace(source, target) if np.random.random() < 0.4 else label
                    for label in labels
                ]

    # If the label is short, add a random suffix
    labels = [
        "%s%s%s"
        % (
            label,
            np.random.choice(["r", "r", "l"]),
            np.random.choice(_SOUNDS["V"]) + "s",
        )
        if len(label) < 5
        else label
        for label in labels
    ]

    # Build the actual labels from genera and epithets and return
    labels = [
        "%s %s" % (_clean_label(genus), _clean_label(epithet).lower())
        for genus, epithet in itertools.zip_longest(labels[:size], labels[size:])
    ]

    return labels
