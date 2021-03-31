import re

from ComputationalBiology.biore.biore_utils import translate


def search(pattern, string, flags=0):
    translated_pattern = translate(pattern)
    return re.search(translated_pattern, string)


def match(pattern, string, flags=0):
    translated_pattern = translate(pattern)
    return re.match(translated_pattern, string)


################

# --------------------------------------------------------------------
# public interface
# TODO: erase all re original documentation
def match(pattern, string, flags=0):
    """Try to apply the pattern at the start of the string, returning
    a Match object, or None if no match was found."""
    translated_pattern = translate(pattern)
    return re.match(translated_pattern, string, flags)


def fullmatch(pattern, string, flags=0):
    """Try to apply the pattern to all of the string, returning
    a Match object, or None if no match was found."""
    translated_pattern = translate(pattern)
    return re.fullmatch(translated_pattern, string, flags)


def search(pattern, string, flags=0):
    """Scan through string looking for a match to the pattern, returning
    a Match object, or None if no match was found."""
    translated_pattern = translate(pattern)
    return re.search(translated_pattern, string, flags)


def sub(pattern, repl, string, count=0, flags=0):
    """Return the string obtained by replacing the leftmost
    non-overlapping occurrences of the pattern in string by the
    replacement repl.  repl can be either a string or a callable;
    if a string, backslash escapes in it are processed.  If it is
    a callable, it's passed the Match object and must return
    a replacement string to be used."""
    translated_pattern = translate(pattern)
    return re.sub(translated_pattern, repl, string, count, flags)


def subn(pattern, repl, string, count=0, flags=0):
    """Return a 2-tuple containing (new_string, number).
    new_string is the string obtained by replacing the leftmost
    non-overlapping occurrences of the pattern in the source
    string by the replacement repl.  number is the number of
    substitutions that were made. repl can be either a string or a
    callable; if a string, backslash escapes in it are processed.
    If it is a callable, it's passed the Match object and must
    return a replacement string to be used."""
    translated_pattern = translate(pattern)
    return re.subn(translated_pattern, repl, string, count, flags)


def split(pattern, string, maxsplit=0, flags=0):
    """Split the source string by the occurrences of the pattern,
    returning a list containing the resulting substrings.  If
    capturing parentheses are used in pattern, then the text of all
    groups in the pattern are also returned as part of the resulting
    list.  If maxsplit is nonzero, at most maxsplit splits occur,
    and the remainder of the string is returned as the final element
    of the list."""
    translated_pattern = translate(pattern)
    return re.split(translated_pattern, string, maxsplit, flags)


def findall(pattern, string, flags=0):
    """Return a list of all non-overlapping matches in the string.
    If one or more capturing groups are present in the pattern, return
    a list of groups; this will be a list of tuples if the pattern
    has more than one group.
    Empty matches are included in the result."""
    translated_pattern = translate(pattern)
    return re.findall(translated_pattern, string, flags)


def finditer(pattern, string, flags=0):
    """Return an iterator over all non-overlapping matches in the
    string.  For each match, the iterator returns a Match object.
    Empty matches are included in the result."""
    translated_pattern = translate(pattern)
    return re.finditer(translated_pattern, string, flags)


def compile(pattern, flags=0):
    "Compile a regular expression pattern, returning a Pattern object."
    translated_pattern = translate(pattern)
    return re.compile(translated_pattern, flags)


def purge():
    "Clear the regular expression caches"
    re.purge()


def template(pattern, flags=0):
    "Compile a template pattern, returning a Pattern object"
    translated_pattern = translate(pattern)
    return re.template(translated_pattern, flags)


def escape(pattern):
    translated_pattern = translate(pattern)
    return re.escape(translated_pattern)
