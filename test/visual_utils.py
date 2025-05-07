# visual_utils.py

from termcolor import colored

class MutationVisualizer:
    @staticmethod
    def highlight_diff(before, after):
        result = []
        for b, a in zip(before, after):
            if b == a:
                result.append(str(a))
            else:
                result.append(colored(str(a), 'red'))
        return ' '.join(result)

    @staticmethod
    def print_encoding_diff(before, after, label="Mutation"):
        print(f"{label} - Original Encoding: {' '.join(map(str, before))}")
        print(f"{label} - Mutated  Encoding: {MutationVisualizer.highlight_diff(before, after)}")
        MutationVisualizer.print_change_summary(before, after)

    @staticmethod
    def print_change_summary(before, after):
        changes = sum(1 for b, a in zip(before, after) if b != a)
        print(f"→ Changed {changes} out of {len(before)} positions.\n")

    @staticmethod
    def print_conformation_pair(conf_before, conf_after, label="Mutation"):
        print(f"\n====== {label} Conformation Before ======")
        conf_before.printAsciiPicture()
        print(f"\n====== {label} Conformation After ======")
        conf_after.printAsciiPicture()
